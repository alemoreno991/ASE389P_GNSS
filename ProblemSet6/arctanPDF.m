classdef arctanPDF
    %ARCTANPDF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mu1_
        mu2_
        sigma1_
        sigma2_
    end
    
    methods
        function obj = arctanPDF(mu1, mu2, sigma1, sigma2)
            obj.mu1_ = mu1;
            obj.mu2_ = mu2;
            obj.sigma1_ = sigma1;
            obj.sigma2_ = sigma2;
        end
        
        function pdf = calculate(obj, mu1, mu2, theta)
            obj.mu1_ = mu1;
            obj.mu2_ = mu2;

            term1 = obj.b(theta) * obj.d(theta) / (sqrt(2*pi)*obj.sigma1_*obj.sigma2_*obj.a(theta)^3);
            term2 = exp(-obj.c()/2)/(pi*obj.sigma1_*obj.sigma2_*obj.a(theta)^2);
            
            Phi1 = normcdf( obj.b(theta)/obj.a(theta));
            Phi2 = normcdf(-obj.b(theta)/obj.a(theta));

            pdf = (term1 * (Phi1 - Phi2) + term2) / cos(theta)^2;
        end

        function a = a(obj, theta)
            a = sqrt( tan(theta)^2/obj.sigma1_^2 + 1/obj.sigma2_^2 );
        end
        
        function b = b(obj, theta)
            b = obj.mu1_/obj.sigma1_^2 * tan(theta) + obj.mu2_/obj.sigma2_^2;
        end
        
        function c = c(obj)
            c = obj.mu1_^2/obj.sigma1_^2 + obj.mu2_^2/obj.sigma2_^2;
        end
        
        function d = d(obj, theta)
            d = exp( (obj.b(theta)^2 - obj.c() * obj.a(theta)^2) / (2* obj.a(theta)^2) );
        end
    end
end

