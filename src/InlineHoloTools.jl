module InlineHoloTools
    using CUDA
    using CUDA.CUFFT
    using Images
    using StatsBase
    using Statistics

    export loadholo, CuTransSqr!, CuTransFunc!, getPhaseRetrievedHolo!
    
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
    function getPhaseRetrievedHolo!(Plane::CuDeviceMatrix{ComplexF32}, img1::CuDeviceMatrix{Float32}, img2::CuDeviceMatrix{Float32}, trans::CuDeviceMatrix{ComplexF32}, transInv::CuDeviceMatrix{ComplexF32}, iterations::Int, datLen::Int)
        compAmp1 = CuArray{ComplexF32}(undef,(datLen,datLen))
        compAmp2 = CuArray{ComplexF32}(undef,(datLen,datLen))
        phi1 = CUDA.ones(datLen,datLen)
        phi2 = CuArray{Float32}(undef,(datLen,datLen))
    
        sqrtImg1 = CuArray{Float32}(undef,(datLen,datLen))
        sqrtImg2 = CuArray{Float32}(undef,(datLen,datLen))
        sqrtImg1 .= sqrt(mean(img1).*255.0)
        sqrtImg2 .= sqrt(mean(img1).*255.0)
        sqrtImg1[div(datLen,4)+1:div(datLen,4)*3,div(datLen,4)+1:div(datLen,4)*3] .= sqrt.(img1.*255.0)
        sqrtImg2[div(datLen,4)+1:div(datLen,4)*3,div(datLen,4)+1:div(datLen,4)*3] .= sqrt.(img2.*255.0)
    
        compAmp1 .= sqrtImg1.*1.0
    
        for itr in 1:iterations
            # STEP 1
            compAmp2 .= CUFFT.ifft(CUFFT.fftshift(CUFFT.fftshift(CUFFT.fft(compAmp1)).*trans))
            phi2 .= angle.(compAmp2) # angle() は複素数の偏角を返す関数。ドットをつけて angle.() で Element-wise な演算をします。
            # STEP 2
            compAmp2 .= sqrtImg2.*exp.(1.0im.*phi2)
            # STEP 3
            compAmp1 .= CUFFT.ifft(CUFFT.fftshift(CUFFT.fftshift(CUFFT.fft(compAmp2)).*transInv))
            phi1 = angle.(compAmp1)
            # STEP 4
            compAmp1 .= sqrtImg1.*exp.(1.0im.*phi1)
        end
    
        Plane .= compAmp1
        return nothing # CUDA 関数は return してはいけません。
    end
end
