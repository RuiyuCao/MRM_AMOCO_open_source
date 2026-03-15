function [g] = LLRp_v3(img, B, lambda, p, useGPU)

[Nx,Ny,Nz,Nt] = size(img);

g = resCalc(img);

    function res = resCalc(img)

        res = img;

        ind = patchInd(img,B);

        nP = numel(ind)/(Nt*B*B*B);

        BLOCK = img(ind);

        BLOCK = reshape(BLOCK,[B*B*B, nP, Nt]);

        BLOCK = permute(BLOCK,[1 3 2]);

        %% -------- SVD shrinkage --------
        if useGPU

            % GPU version
            for ip = 1:nP
                BLOCK(:,:,ip) = svt(BLOCK(:,:,ip), lambda, p);
            end

        else

            % CPU parallel version
            parfor ip = 1:nP
                BLOCK(:,:,ip) = svt(BLOCK(:,:,ip), lambda, p);
            end

        end

        BLOCK = permute(BLOCK,[1 3 2]);

        tempres = res;

        tempres(ind(:,1)) = BLOCK;

        res = tempres;

    end

end
