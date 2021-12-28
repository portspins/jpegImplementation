clear all;

testImgOne = imread("Project_Code_Images/alu.tif");
compressedImgOne = jpeg(testImgOne);
imshow(compressedImgOne);
imwrite(compressedImgOne, "output_alu_8x8.png");
psnrOne = calculatePSNR(testImgOne, compressedImgOne);
disp(psnrOne);
calculatePWE(testImgOne, compressedImgOne);

testImgTwo = imread("Project_Code_Images/tulips.png");
compressedImgTwo = jpeg(testImgTwo);
imshow(compressedImgTwo);
imwrite(compressedImgTwo, "output_tulips_8x8.png");
psnrTwo = calculatePSNR(testImgTwo, compressedImgTwo);
disp(psnrTwo);

% jpegTest(testImgOne, "alu");
% jpegTest(testImgTwo, "tulip");

% Applies JPEG process to RGB image
% Param: RGB image
% Return: Compressed RGB Image
function compressedImg = jpeg(inputImg)
    %initialize quantization matrices
    lQuant = [16 11 10 16 24 40 51 61;
              12 12 14 19 26 58 60 55;
              14 13 16 24 40 57 69 56;
              14 17 22 29 51 87 80 62;
              18 22 37 56 68 109 103 77;
              24 35 55 64 81 104 113 92;
              49 64 78 87 103 121 120 101;
              72 92 95 98 112 100 103 99];
    
    cQuant = [17 18 24 47 99 99 99 99;
              18 21 26 66 99 99 99 99;
              24 26 56 99 99 99 99 99;
              47 66 99 99 99 99 99 99;
              99 99 99 99 99 99 99 99;
              99 99 99 99 99 99 99 99;
              99 99 99 99 99 99 99 99;
              99 99 99 99 99 99 99 99];

    %convert to ycbcr and apply chroma 4:2:0 subsampling
    Cb = double(convertRGBToYCbCr(inputImg));
    Cb = chroma420(Cb);
    
    %variable for max row and column sizes (width and height of image)
    rowsize = size(Cb, 1);
    colsize = size(Cb, 2);
    i = 0;
    j = 0;
    %find out if image makes 8x8 cleanly, or if there's remainder
    if (mod(rowsize, 8) ~= 0)
        i = 8 - mod(rowsize, 8);
    end
    if (mod(colsize, 8) ~= 0)
        j = 8 - mod(colsize, 8);
    end
    
    %check i and j, if they are not 0 we need a new (bigger) matrix
    if (i ~= 0 || j ~= 0)
        newrowsize = rowsize+i;
        newcolsize = colsize+j;
        newcb = Cb;
        newcb(rowsize+1:newrowsize, colsize+1:newcolsize, 1:3) = 0;
    else
        %if matrix is perfectly divisible by 8x8 squares
        newrowsize = rowsize;
        newcolsize = colsize;
        newcb = Cb;
    end
    
    %loop through image, in increments of 8
    for  x = 1 : 8 : newrowsize
        for  y = 1 : 8 : newcolsize
            for z = 1 : 3
                sub = newcb(x:x+7, y:y+7, z);
                subDct = dct(sub);
                quant = cQuant;
                if (z == 1)
                    quant = lQuant;
                end
                quantized = quantizeDct(subDct, quant);
                dequantized = dequantizeDct(quantized, quant);
                iDct = idct(dequantized);
                newcb(x:x+7, y:y+7, z) = iDct;
            end
        end
    end
    Cb = newcb(1:rowsize, 1:colsize, 1:3);
    compressedImg = convertYCbCrToRGB(Cb);
end

% Converts RGB image to YCbCr image
% Param: An RGB image
% Return: A YCbCr image
function convertedImg = convertRGBToYCbCr(rgbImg)
    rgbImg = double(rgbImg) / 255; % Scale the RGB image
    convertedImg = rgbImg; % Initialize size of new image
    convertMatrix = [ 0.299 0.587 0.114; -0.168736 -0.331264 0.5; 0.5 -0.418688 -0.081312 ]; % Initialize conversion matrix
    [rows, cols, ~] = size(rgbImg); % Get size of image
    
    for x = 1:rows
        for y = 1:cols
            % Apply the conversion formula to each pixel
            convertedImg(x, y, 1:3) = convertMatrix * [ rgbImg(x, y, 1); rgbImg(x, y, 2); rgbImg(x, y, 3) ] + [ 0; 0.5; 0.5 ];
        end
    end
    
    convertedImg = uint8(convertedImg * 255);
end

% Converts YCbCr image to RGB image
% Param: A YCbCr image
% Return: An RGB image
function convertedImg = convertYCbCrToRGB(ycbcrImg)
    ycbcrImg = double(ycbcrImg) / 255; % Scale the image
    convertedImg = ycbcrImg; % Initialize size of new image
    convertMatrix = [ 0.299 0.587 0.114; -0.168736 -0.331264 0.5; 0.5 -0.418688 -0.081312 ]; % Initialize conversion matrix
    [rows, cols, ~] = size(ycbcrImg); % Get size of image
    
    for x = 1:rows
        for y = 1:cols
            % Apply the conversion formula to each pixel
            convertedImg(x, y, 1:3) = convertMatrix \ ([ ycbcrImg(x, y, 1); ycbcrImg(x, y, 2); ycbcrImg(x, y, 3) ] - [ 0; 0.5; 0.5 ]);
        end
    end
    
    convertedImg = uint8(convertedImg * 255);
end

% Applies Chroma 4:2:0 subsampling to YCbCr image
% Param: A YCbCr image
% Return: Subsampled YCbCr Image
function subsampled2 = chroma420(ycbcrImg)
    subsampled2 = ycbcrImg; % Initialize output image
    [rows, cols, ~] = size(ycbcrImg); % Get image size
    
    for x = 1:2:rows % For every other row in image
        for y = 1:2:cols % And every other pixel in those rows
            % Update so each 2 x 2 block of pixels has same CbCr data
            cbVal = ycbcrImg(x, y, 2);
            crVal = ycbcrImg(x, y, 3);
            if (x < rows && y < cols)
                subsampled2(x:(x + 1), y:(y + 1), 2) = cbVal;
                subsampled2(x:(x + 1), y:(y + 1), 3) = crVal;
            % Handle special case of last row
            elseif (x == rows && y < cols)
                subsampled2(x, y:(y + 1), 2) = cbVal;
                subsampled2(x, y:(y + 1), 3) = crVal;
            % Handle special case of last column
            elseif (x < rows && y == cols)
                subsampled2(x:(x + 1), y, 2) = cbVal;
                subsampled2(x:(x + 1), y, 3) = crVal;
            end
        end
    end
    
end

% Applies 2D DCT to YCbCr image
% Param: An 8x8 Y, Cb, or Cr matrix
% Return: The image with DCT applied
function dctImg = dct(ycbcrImg)
    [rows, cols] = size(ycbcrImg);
    if (rows ~= 8 || cols ~= 8)
        return;
    end
    dctImg = zeros(8,8);
    for u = 0:1:7
        uVal = C(u);
        for v = 0:1:7
            vVal = C(v);
            dctVal = 0;
            for i = 0:1:7
                iVal = cos((2*i+1)*u*pi/16);
                for j = 0:1:7
                    jVal = cos((2*j+1)*v*pi/16);
                    dctVal = dctVal + iVal * jVal * ycbcrImg(i+1,j+1);
                end
            end
            dctVal = dctVal * (uVal * vVal / 4);
            dctImg(u+1,v+1) = round(dctVal);
        end
    end
end

% Applies 2D inverse DCT to YCbCr image
% Param: An 8x8 YCbCr image with DCT applied
% Return: The image with inverse DCT applied
function idctImg = idct(dctImg)
    [rows, cols, ~] = size(dctImg);
    if (rows ~= 8 || cols ~= 8)
        return;
    end
    idctImg = zeros(8,8);
    for i = 0:1:7
        for j = 0:1:7
            idctVal = 0;
            for u = 0:1:7
                uVal = C(u);
                iVal = cos((2*i+1)*u*pi/16);
                for v = 0:1:7
                    vVal = C(v);
                    jVal = cos((2*j+1)*v*pi/16);
                    idctVal = idctVal + iVal * jVal * dctImg(u+1,v+1) * (uVal * vVal / 4);
                end
            end
            idctImg(i+1,j+1) = idctVal;
        end
    end
end

% Helper function for DCT
% Param: Index
% Return: Corresponding C value
function cVal = C(u)
    if (u == 0)
        cVal = 1 / sqrt(2);
    else
        cVal = 1;
    end
end

% Quantizes DCT block
% Params: An 8x8 DCT block, the 8x8 quantization matrix
% Return: The quantized block
function quantized = quantizeDct(dctBlock, qBlock)
    quantized = round(dctBlock ./ qBlock);
end

% Dequantizes quantized DCT block
% Params: A quantized 8x8 DCT block, the 8x8 quantization matrix
% Return: The dequantized block
function dequantized = dequantizeDct(dctBlock, qBlock)
    dequantized = dctBlock .* qBlock;
end

% Calculates the peak signal-to-noise ratio
% Params: An original reference image, a compressed image
% Return: The peak signal-to-noise ratio
function psnr = calculatePSNR(originalImg, compressedImg)
    % Retrieve width and height of images
    rows = size(originalImg, 1);
    cols = size(originalImg, 2);

    % Calculate mean-square error of each color channel
    mseR = sum(sum((double(originalImg(:,:,1)) - double(compressedImg(:,:,1))) .^ 2)) / (rows * cols);
    mseG = sum(sum((double(originalImg(:,:,2)) - double(compressedImg(:,:,2))) .^ 2)) / (rows * cols);
    mseB = sum(sum((double(originalImg(:,:,3)) - double(compressedImg(:,:,3))) .^ 2)) / (rows * cols);
    
    % Calculate mean-square error for all channels by averaging
    mse = mean([mseR, mseG, mseB]);
    
    % Calculate peak signal-to-noise ratio
    psnr = 20 * log10( 255 / sqrt(mse));
end

% Calculates pixelwise error
% Params: An original reference image, a compressed image
function pwe = calculatePWE(originalImg, compressedImg)
    error = compressedImg - originalImg;
    imagesc(error);
end

% Applies JPEG process to first 8x8 block of RGB image
% Param: RGB image
function compressedImg = jpegTest(inputImg, str)
    %initialize quantization matrices
    lQuant = [16 11 10 16 24 40 51 61;
              12 12 14 19 26 58 60 55;
              14 13 16 24 40 57 69 56;
              14 17 22 29 51 87 80 62;
              18 22 37 56 68 109 103 77;
              24 35 55 64 81 104 113 92;
              49 64 78 87 103 121 120 101;
              72 92 95 98 112 100 103 99];
    
    newImg = inputImg(1:8, 1:8, 1:3);
    imwrite(newImg, append(str, '_', "RGBbefore.png"));
    %convert to ycbcr and apply chroma 4:2:0 subsampling
    Cb = double(convertRGBToYCbCr(newImg));
    
    disp(Cb(1:8, 1:8, 2:3));
    Cb = chroma420(Cb);

    sub = Cb(1:8, 1:8, 1);
    disp(sub);
    subDct = dct(sub);
    disp(subDct);
    quantized = quantizeDct(subDct, lQuant);
    disp(quantized);
    dequantized = dequantizeDct(quantized, lQuant);
    disp(dequantized);
    iDct = idct(dequantized);
    disp(round(iDct));
end
