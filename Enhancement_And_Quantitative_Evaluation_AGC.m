%This code enhances images using AGC algorithm and computes rms and entropy
%of the enhanced images

clear;
clc;

path1=strcat(pwd,'\Database');
listing=dir(path1);
SZ=size(listing);
counter=0;

rms_I_AGC=zeros(1,SZ(1)-2); %Because the first two are always non-images
ent_I_AGC=zeros(1,SZ(1)-2);

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
        if(NN==3) %RGB images
            B=rgb2hsv(Im);
            A=B(:,:,3);
        else
            A=double(Im);
            A=A/255;
        end

        M=mean(A(:));
        SD=std(A(:));

        %Computation of gamma
        diff=(M+(2*SD))-(M-(2*SD));
        if(diff<=0.333333)
            gm=log2(SD);
            gm=-gm;
        else
            gm=exp((1-(M+SD))/2);
        end
    
    
        SS1=size(A);
        Iout=zeros(SS1(1),SS1(2));
    
        trm=(0.5-M);
        if(trm<0)
            trm2=0;
        else
            trm2=1;
        end
        sm=0;
        for i=1:SS1(1)
            for j=1:SS1(2)
                %Computation of c
                k=(A(i,j)^gm)+(1-(A(i,j))^gm)*(M^gm);           
                CC=1/(1+(trm2*(k-1)));
                Iout(i,j)=CC*A(i,j)^gm;        
            end
        end
    
        if(NN==3) %For RGB images
            ImOut1=zeros(SS(1),SS(2),SS(3));
            ImOut1(:,:,1:2)=B(:,:,1:2);    
            ImOut1(:,:,3)=Iout;
            ImOut1=hsv2rgb(ImOut1);
            ImOut1=255*ImOut1;
            ImOut1=uint8(ImOut1);
        else %Non-RGB images       
            ImOut1=uint8(255*(Iout));
        end

        %Uncomment the below three lines to view the original and enhanced images
        %Since the database is huge, uncommenting is suggested only for
        %selected images (step running).

        %imshow(Im);
        %figure;
        %imshow(ImOut1);
    
        %rms computation                
        Im_I=double(ImOut1);               
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
            rms_I_AGC(counter)=(SM_I/(SS(1)*SS(2)*SS(3)))^0.5;
        else %Non-RGB images
            for i=1:SS(1)
                for j=1:SS(2)                 
                    dif_I=(M_I-Im_I(i,j))^2;
                    SM_I=SM_I+dif_I;                    
                end
            end
            rms_I_AGC(counter)=(SM_I/(SS(1)*SS(2)))^0.5;
        end      
    
        %Entropy computation
        if(NN==3)%RGB images
            ent_I_R=entropy(ImOut1(:,:,1));
            ent_I_G=entropy(ImOut1(:,:,2));
            ent_I_B=entropy(ImOut1(:,:,3));
            ent_I_AGC(counter)=(ent_I_R+ent_I_G+ent_I_B)/3;
        else
            ent_I_AGC(counter)=entropy(ImOut1);
        end 
    end
end

%Mean entropy and rms of the database
ent_mean_AGC=mean(ent_I_AGC);
rms_mean_AGC=mean(rms_I_AGC);

save Mean_Entropy_AGC ent_mean_AGC;
save Mean_RMS_AGC rms_mean_AGC;