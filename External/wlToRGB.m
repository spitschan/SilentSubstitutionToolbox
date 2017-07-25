function RGB = Wavelength_to_RGB(wl)

if nargin<1
    wl = 650;
else
    if (wl<380 || wl>780)
        error('ERROR! (Given wavelength out of visible range [380 780])')
        return
    end
end

M  = 780-380;
CV = zeros(1,M,3);

for ii=1:M

    % WAVELENGTH = WL
    WL = 380 + real(ii * 400 / M);

    if (WL>=380 && WL<=440)
        R = -1.*(WL-440)/(440-380);
        G = 0;
        B = 1;
    end
    if (WL>=440 && WL<=490)
        R = 0;
        G = (WL-440)/(490-440);
        B = 1;
    end
    if (WL>=490 && WL<=510)
        R = 0;
        G = 1;
        B = -1.*(WL-510)/(510-490);
    end
    if (WL>=510 && WL<=580)
        R = (WL-510)/(580-510);
        G = 1;
        B = 0;
    end
    if (WL>=580 && WL<=645)
        R = 1;
        G = -1.*(WL-645)/(645-580);
        B = 0;
    end
    if (WL>=645 && WL<=780)
        R = 1;
        G = 0;
        B = 0;
    end

    % LET THE INTENSITY SSS FALL OFF NEAR THE VISION LIMITS
    if (WL>700)
        SSS=0.3+0.7* (780-WL)/(780-700);
    elseif (WL<420)
        SSS=.3+.7*(WL-380)/(420-380);
    else
        SSS=1;
    end

    CV(1,ii+380,1) = SSS*R;
    CV(1,ii+380,2) = SSS*G;
    CV(1,ii+380,3) = SSS*B;
end


% CALCULATE R-G-B VALUE OF THE GIVEN WAVELENGTH
RGB = round([CV(1,round(wl),1) CV(1,round(wl),2) CV(1,round(wl),3)].*255); 