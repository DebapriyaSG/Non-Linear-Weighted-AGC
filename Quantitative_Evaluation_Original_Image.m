%This code computes rms and entropy of images without any enhancement

clear;
clc;

path1=strcat(pwd,'\Database');
listing=dir(path1);
SZ=size(listing);
counter=0;

rms_I_Original=zeros(1,SZ(1)-2); %Because the first two are always non-images
ent_I_Original=zeros(1,SZ(1)-2);

for ii=1:SZ(1)  
    f=listing(ii).name;
    byt=listing(ii).bytes;         

    if(byt~=0)
        counter=counter+1;
        file=strcat(path1,'\',f);     
       
        Im=imread(file);
        %Uncomment the following line for viewing the images. Since the database is huge, uncommenting is suggested only for
        %selected images (step running).
        %imshow(Im);

        if(islogical(Im))
            Im=uint8(255*Im);
        end
        SS=size(Im);
        NN=ndims(Im);
        
       %rms computation    
        Im_I=double(Im);               
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
            rms_I_Original(counter)=(SM_I/(SS(1)*SS(2)*SS(3)))^0.5;
        else %Non-RGB images
            for i=1:SS(1)
                for j=1:SS(2)                 
                    dif_I=(M_I-Im_I(i,j))^2;
                    SM_I=SM_I+dif_I;                    
                end
            end
            rms_I_Original(counter)=(SM_I/(SS(1)*SS(2)))^0.5;
        end      
    
        if(NN==3) %RGB images
            ent_I_R=entropy(Im(:,:,1));
            ent_I_G=entropy(Im(:,:,2));
            ent_I_B=entropy(Im(:,:,3));
            ent_I_Original(counter)=(ent_I_R+ent_I_G+ent_I_B)/3;
        else %Non-RGB images
            ent_I_Original(counter)=entropy(Im);
        end 
            
     end         
 end   
    
    %Mean entropy and rms of the database
    ent_mean_original=mean(ent_I_Original);
    rms_mean_original=mean(rms_I_Original);   
    
    
    save Mean_Entropy_Original ent_mean_original;
    save Mean_RMS_Original rms_mean_original;