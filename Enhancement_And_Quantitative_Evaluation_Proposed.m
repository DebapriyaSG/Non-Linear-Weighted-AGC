%This code enhances images using the proposed non-linear weighted AGC algorithm
%and computes rms and entropy of the enhanced images

clear;
clc;

path1=strcat(pwd,'\Database');
listing=dir(path1);
SZ=size(listing);
counter=0;

rms_I=zeros(1,SZ(1)-2); %Because the first two are always non-images
ent_I=zeros(1,SZ(1)-2);

for ii=1:SZ(1)  
    f=listing(ii).name;
    byt=listing(ii).bytes;         

    if(byt~=0)
        counter=counter+1;
        file=strcat(path1,'\',f);

        Im=imread(file);
        if(islogical(Im))
            Im=uint8(255*Im);
        end
        SS=size(Im);
        NN=ndims(Im);
        if(NN==3)%RGB images
            B=rgb2hsv(Im);
            A=B(:,:,3);
        else
            A=double(Im);
            A=A/255;
        end

        M=mean(A(:));
        SD=std(A(:));

        %Computation of gamma_nu
        pr_g=2.47;
        gm=((2*SD)^(-log2((2*SD)^pr_g)))*(exp((1-(M+SD))/2))+((1-((2*SD)^(-log2((2*SD)^pr_g))))*((-1)*log2(SD)));

        SS1=size(A);
        Iout=zeros(SS1(1),SS1(2));
        pr_c=9.96;

        for i=1:SS1(1)
            for j=1:SS1(2)

                %Computation of c_nu
                k=M^(-log2(M^pr_c))+((1-(M^(-log2(M^pr_c))))*((A(i,j)^gm)+(1-(A(i,j))^gm)*(M^gm)));
                Iout(i,j)=(A(i,j)^gm)/k;
            end
        end

        if(NN==3) %For RGB images
            ImOut=zeros(SS(1),SS(2),SS(3));
            ImOut(:,:,1:2)=B(:,:,1:2);    
            ImOut(:,:,3)=Iout;
            ImOut=hsv2rgb(ImOut);
            ImOut=255*ImOut;
            ImOut=uint8(ImOut);
        else %Non-RGB images           
            ImOut=uint8(255*(Iout));
        end
        
        %Uncomment the below four lines to view the original and enhanced images
        %Since the database is huge, uncommenting is suggested only for
        %selected images (step running).
        
        %figure;
        %imshow(Im);
        %figure;
        %imshow(ImOut);

        %rms computation  
        Im_I=double(ImOut);               
        Im_I=Im_I/255;                
        M_I=mean(Im_I(:));                
        SM_I=0;
        if(NN==3) %RGB images
            for i=1:SS(1)
                for j=1:SS(2)
                    for kk=1:SS(3)                                
                        dif_I=(M_I-Im_I(i,j,kk))^2;                                
                        SM_I=SM_I+dif_I;
                    end
                end
            end
            rms_I(counter)=(SM_I/(SS(1)*SS(2)*SS(3)))^0.5;
        else %Non-RGB images
            for i=1:SS(1)
                for j=1:SS(2)                 
                    dif_I=(M_I-Im_I(i,j))^2;
                    SM_I=SM_I+dif_I;                    
                end
            end
            rms_I(counter)=(SM_I/(SS(1)*SS(2)))^0.5;
        end         

        
        %Entropy computation

        if(NN==3) %RGB images 
            ent_I_R=entropy(ImOut(:,:,1));
            ent_I_G=entropy(ImOut(:,:,2));
            ent_I_B=entropy(ImOut(:,:,3));
            ent_I(counter)=(ent_I_R+ent_I_G+ent_I_B)/3;
        else %Non-RGB images
            ent_I(counter)=entropy(ImOut);
        end 
    end
end

%Mean entropy and rms of the database
ent_mean_proposed=mean(ent_I);
rms_mean_proposed=mean(rms_I);

save Mean_Entropy_Proposed ent_mean_proposed;
save Mean_RMS_Proposed rms_mean_proposed;
        
       

