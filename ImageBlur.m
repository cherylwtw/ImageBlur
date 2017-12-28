function main()
    clear all;
    close all;
    
    sd = 4;
    gaussianKernel = GetTwoDGaussian(sd);
    % for drawing
    steps = sd*2*2+1;
    xValues = linspace(-2*sd,2*sd,steps);
    yValues = linspace(-2*sd,2*sd,steps);
    
    % draw Gaussian
    figure;
    surf(xValues,yValues,gaussianKernel);
    
    gaussianBlurredImage = GaussianBlur(sd);
    
    fourierTransformBlurred = FourierTransformBlur(sd);
    fourierTransformBlurredUint8 = uint8(fourierTransformBlurred);
    croppedFourierTransformBlurred = imcrop(fourierTransformBlurredUint8, [9, 9, 426, 426]);
    
    % show result from both blurs
    %figure;
    %imshow(gaussianBlurredImage);
    %title("Padded Gaussian Blurred Image");
    
    croppedGaussinImage = imcrop(gaussianBlurredImage, [9, 9, 426, 426]);
    figure;
    imshow(croppedGaussinImage);
    title("Gaussian Blurred Image");
    
    %figure;
    %imshow(fourierTransformBlurredUint8);
    %title("Padded Fourier Transformed Blurred Image");
    
    figure;
    imshow(uint8(croppedFourierTransformBlurred));
    title("Fourier Transformed Blurred Image");
    mseBetweenGaussionBlurAndFourierTransformBlur = immse (croppedGaussinImage, croppedFourierTransformBlurred)
    
    stepCount = (sd^2)/2;
    heatBlurred = HeatFunctionBlur(stepCount);
    heatBlurredUint8 = uint8(heatBlurred);
    figure;
    imshow(heatBlurredUint8);
    title("Heat Equation Blurred Image");
    mseBetweenGaussianBlurAndHeatEquationBlur = immse(croppedGaussinImage,heatBlurredUint8)
    %immse(croppedFourierTransformBlurred,heatBlurredUint8)
end

% Draw 2D Gaussian
function gaussianKernel = GetTwoDGaussian(sd)
    steps = sd*2*2+1;
    gaussianKernel = zeros(steps, steps);

    xValues = linspace(-2*sd,2*sd,steps);
    yValues = linspace(-2*sd,2*sd,steps);

    for xIdx = 1:length(xValues)
        for yIdx = 1:length(yValues)
            x = xValues(xIdx);
            y = xValues(yIdx);
            r = double(sqrt(x.^2 + y.^2));
            g = (1/(sqrt(2*pi).*sd.^2))*exp(-(r.^2)/(2.*(sd.^2)));
            gaussianKernel(xIdx, yIdx) = g;
        end
    end
    
   % devide by the sum, so when it convolutes with the image, the pixel
   % value does not exceed the maximum
   gaussianKernel=gaussianKernel./sum(gaussianKernel(:)); 
end

% Gaussian Blur
function gaussianBlurredImage = GaussianBlur(sd)
    kernelSize = sd*2*2+1;

    % build gaussian kernel
    gaussianKernel = GetTwoDGaussian(sd);

    % get grayscale picture
    mitPicFileName = strcat('mit','.jpg');
    mitColorImage = imread(mitPicFileName);
    mitGreyImage = rgb2gray(mitColorImage);
    
    % add zero paddings
    paddedGreyImage = padarray(mitGreyImage, [8 8], 0, 'both');
    [greyRows, greyColumns, numberOfGreyIntensityChannels] = size(paddedGreyImage);

    % create a new image to store new blurred pixel value
    blurredGreyImage=uint8(zeros(greyRows,greyColumns,1));

    for row = 1:greyRows
        for col = 1:greyColumns
            kernelRowStart = 1;
            kernelRowEnd = kernelSize;
            kernelColStart = 1;
            kernelColEnd = kernelSize;

            % only assign blurredGreyImage with the ones that have kernel
            % fully covered
            if ( row < greyRows - kernelSize+2 && col < greyColumns - kernelSize+2)
                
            intensityTotal = 0;
            for kRow = kernelRowStart:kernelRowEnd
                for kCol = kernelColStart:kernelColEnd
                    intensityTotal = intensityTotal + ...
                        double(paddedGreyImage(row+kRow-1,...
                        col+kCol-1, 1))*...
                        gaussianKernel(kRow,kCol);
                end
            end
            blurredPixel = double(intensityTotal);
            % assign the center pixel of where the kernel covers
            blurredGreyImage(row+sd*2,col+sd*2,1) = uint8(blurredPixel);
            end 
        end
    end
    gaussianBlurredImage = blurredGreyImage;
end

function fourierTransformBlurred = FourierTransformBlur(sd)
    % fourier transform of image
    mitPicFileName = strcat('mit','.jpg');
    mitColorImage = imread(mitPicFileName);
    mitGreyImage = rgb2gray(mitColorImage);
    
    % 8 padding on each side like the previous method
    paddedSiz = 427+8;
    [imageRows, imageCols, layers] = size(mitGreyImage);
    imageTransformed = fft2(double(mitGreyImage),paddedSiz, paddedSiz);
    
    % fourier tranform of gaussian
    gaussianKernel = GetTwoDGaussian(sd);
    gaussianTransformed = fft2(gaussianKernel, paddedSiz, paddedSiz);

    % scaler multiply of two fourier transformed
    fourierTransformed=double(imageTransformed.*gaussianTransformed);
    
    % get the inverse fourier transform of the scaler multiply result
    inversedFourierTransformed = ifft2(double(fourierTransformed));
    fourierTransformBlurred = inversedFourierTransformed;
end

function heatBlurredImage = HeatFunctionBlur(step)
    stepSize = 1;
    
    mitPicFileName = strcat('mit','.jpg');
    mitColorImage = imread(mitPicFileName);
    originalImage = rgb2gray(mitColorImage);
    originalImage = double(originalImage);
    
    [greyRows, greyColumns, numberOfGreyIntensityChannels] = size(originalImage);
    
    imageAfterStep =double(zeros(greyRows,greyColumns));
    
    for s = 1: step
        [ix, iy] = gradient(originalImage);
        [ixx,h] = gradient(ix);
        [h,iyy] = gradient(iy);
        
        g = ixx+iyy;
        
        imageAfterStep = stepSize*g+double(originalImage);
        
        %imshow(mat2gray(originalImage));
        originalImage = imageAfterStep;
    end
    
    heatBlurredImage = imageAfterStep;
end
    
