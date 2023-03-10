module InlineHoloTools
    using CUDA
    using CUDA.CUFFT
    using Images
    using StatsBase
    using Statistics

    export loadholo, CuTransSqr!, CuTransFunc!, getPhaseRetrievedHolo!, getReconstvolfromHolo!
    
    """
    Load a grayscale hologram and return as a Matrix
    """
    function loadholo(path::String)
        out = Float32.(channelview(Gray.(load(path))))
    end

    """
    **This is a CUDA kernel**\n
    Calculate of the root sign matrix for diffraction integrals by Fourier transform and transfer function (Kreis, et al., 1997). The center of the even-numbered array, i.e., `div(datLen,2)+1` th element, is `(1,1)` from the coordinate origin.
    """
    function CuTransSqr!(Plane::CuDeviceMatrix{Float32}, datLen::Int, wavLen::Float32, dx::Float32)
        x = (blockIdx().x-1)*blockDim().x + threadIdx().x
        y = (blockIdx().y-1)*blockDim().y + threadIdx().y
        if x <= datLen && y <= datLen
            Plane[x,y] = 1.0 - ((x-datLen/2)*wavLen/datLen/dx)^2 - ((y-datLen/2)*wavLen/datLen/dx)^2
        end
        return
    end

    """
    **This is a CUDA kernel**\n
    Define and calculate the transfer function. `d_sqrPart` needs to be calculated with `CuTransSqr!()` in advance.
    """
    function CuTransFunc!(Plane::CuDeviceMatrix{ComplexF32}, z0::Float32, wavLen::Float32, datLen::Int, d_sqrPart::CuDeviceMatrix{Float32})
        x = (blockIdx().x-1)*blockDim().x + threadIdx().x
        y = (blockIdx().y-1)*blockDim().y + threadIdx().y
        if x <= datLen && y <= datLen
            Plane[x,y] = exp(2im*pi*z0/wavLen*sqrt(d_sqrPart[x,y]))
        end
        return
    end

    """
    **This is *NOT* a CUDA kernel**\n
    Retrieve the phase information of recorded light with two holograms `img1` and `img2` by using Gerchberg-Saxton algorithm. The phase-retrieved plane is given as `Plane`. Forward and reverse transfer functions (`trans` and `transInv`) based on the distance between the two recorded holograms are needed, respectively. `iterations` specifies the number of times to perform iterative phase recovery based on the GS algorithm.
    """
    function getPhaseRetrievedHolo!(Plane::CuArray{ComplexF32, 2, CUDA.Mem.DeviceBuffer}, img1::CuArray{Float32, 2, CUDA.Mem.DeviceBuffer}, img2::CuArray{Float32, 2, CUDA.Mem.DeviceBuffer}, trans::CuArray{ComplexF32, 2, CUDA.Mem.DeviceBuffer}, transInv::CuArray{ComplexF32, 2, CUDA.Mem.DeviceBuffer}, iterations::Int, imgLen::Int, padFlag::Bool)
        datLen::Int = 0
        if padFlag
            datLen = imgLen*2
        else
            datLen = imgLen
        end

        compAmp1 = CuArray{ComplexF32}(undef,(datLen,datLen))
        compAmp2 = CuArray{ComplexF32}(undef,(datLen,datLen))
        phi1 = CUDA.ones(datLen,datLen)
        phi2 = CuArray{Float32}(undef,(datLen,datLen))
    
        sqrtImg1 = CuArray{Float32}(undef,(datLen,datLen))
        sqrtImg2 = CuArray{Float32}(undef,(datLen,datLen))
        sqrtImg1 .= sqrt(mean(img1).*255.0)
        sqrtImg2 .= sqrt(mean(img1).*255.0)

        if padFlag
            sqrtImg1[div(datLen,4)+1:div(datLen,4)*3,div(datLen,4)+1:div(datLen,4)*3] .= sqrt.(img1.*255.0)
            sqrtImg2[div(datLen,4)+1:div(datLen,4)*3,div(datLen,4)+1:div(datLen,4)*3] .= sqrt.(img2.*255.0)
        else
            sqrtImg1 .= sqrt.(img1.*255.0)
            sqrtImg2 .= sqrt.(img2.*255.0)
        end
    
        compAmp1 .= sqrtImg1.*1.0
    
        for itr in 1:iterations
            # STEP 1
            compAmp2 .= CUFFT.ifft(CUFFT.fftshift(CUFFT.fftshift(CUFFT.fft(compAmp1)).*trans))
            phi2 .= angle.(compAmp2) # angle() ???????????????????????????????????????????????????????????? angle.() ??? Element-wise ????????????????????????
            # STEP 2
            compAmp2 .= sqrtImg2.*exp.(1.0im.*phi2)
            # STEP 3
            compAmp1 .= CUFFT.ifft(CUFFT.fftshift(CUFFT.fftshift(CUFFT.fft(compAmp2)).*transInv))
            phi1 = angle.(compAmp1)
            # STEP 4
            compAmp1 .= sqrtImg1.*exp.(1.0im.*phi1)
        end
    
        Plane .= compAmp1
        return nothing # CUDA ????????? return ???????????????????????????
    end

    """
    **This is *NOT* a CUDA kernel**\n
    The reconstructed volume is obtained from the complex amplitude. The complex amplitude is obtained by phase reconstruction using the GS method or by taking the square root of the intensity of the hologram image and casting the array in `ComplexF32`. A total of `recItr` slices can be obtained for each dz from the plane defined by `frontz`.
    """
    function getReconstvolfromHolo!(vol::CuArray{Float32, 3, CUDA.Mem.DeviceBuffer}, holo::CuArray{ComplexF32, 2, CUDA.Mem.DeviceBuffer}, transfront::CuArray{ComplexF32, 2, CUDA.Mem.DeviceBuffer}, transdz::CuArray{ComplexF32, 2, CUDA.Mem.DeviceBuffer}, recItr::Int, datLen::Int)
        holotmp = CuArray{ComplexF32}(undef,datLen,datLen)
        holo .= CUFFT.fftshift(CUFFT.fft(holo)) .* transfront
        holotmp .= CUFFT.ifft(CUFFT.fftshift(holo))
        vol[1,:,:] .= Float32.(abs.(holotmp .* conj.(holotmp)))

        for itr in 2:recItr
            holo .= holo .* transdz
            holotmp .= CUFFT.ifft(CUFFT.fftshift(holo))
            vol[itr,:,:] .= Float32.(abs.(holotmp .* conj.(holotmp)))
        end
    end
    return nothing
end
