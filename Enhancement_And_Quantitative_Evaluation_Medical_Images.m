%This code enhances MR and PET slices of brain using non-linear weighted
%AGC algorithm, stores the enhanced images, and compares entropy and rms of
%enhanced and original images

clear;
clc;

path1=strcat(pwd,'\Medical_Data');
listing=dir(path1);
SZ=size(listing);

rms_proposed=zeros(1,SZ(1)-2); %Because the first two are always non-images
ent_proposed=zeros(1,SZ(1)-2);
rms_original=zeros(1,SZ(1)-2);
ent_original=zeros(1,SZ(1)-2);
counter=0;

path2=strcat(pwd,'\Output_Medical_Data');

if ~exist(path2, 'dir')
       mkdir(path2)
end


for ii=1:SZ(1)  
    f=listing(ii).name;
    byt=listing(ii).bytes;         

    if(byt~=0)        
        counter=counter+1;
        file=strcat(path1,'\',f);
        
        A=dicomread(file);        
        B=A;        
        A=double(A);
        l=min(A(:));
        h=max(A(:));
        SS=size(A);

        A=A-l;        
        A=A/h;

        M=mean(A(:));
        SD=std(A(:));  
        

       % Original image rms and entropy for comparison
        SM_I=0;
        for i=1:SS(1)
            for j=1:SS(2)                 
                dif_I=(M-A(i,j))^2;
                SM_I=SM_I+dif_I;                    
            end
        end
        rms_original(counter)=(SM_I/(SS(1)*SS(2)))^0.5;
        ent_original(counter)=entropy(uint16(B)); %Entropy
       
        %Enhancement with proposed algorithm   
    
        %Computation of gamma_nu
        gm=((2*SD)^(-log2((2*SD)^2.47)))*(exp((1-(M+SD))/2))+((1-((2*SD)^(-log2((2*SD)^2.47))))*((-1)*log2(SD)));
        Iout=zeros(SS(1),SS(2));   
    
        for i=1:SS(1)
            for j=1:SS(2)  
                %Computation of c_nu
                k=M^(-log2(M^9.96))+((1-(M^(-log2(M^9.96))))*((A(i,j)^gm)+(1-(A(i,j))^gm)*(M^gm)));
                Iout(i,j)=(A(i,j)^gm)/k;        
            end
        end
        
        ImOut=int16(h*(Iout));
    
        file_name=strcat('Enhanced_',num2str(counter),'.dcm'); 
        path3=strcat(path2,'\',file_name);      
        dicomwrite(ImOut,path3);     
        
        %rms computation        
        Im_I=double(ImOut);        
        Im_I=Im_I/h;                
        M_I=mean(Im_I(:));                
        SM_I=0;
    
        for i=1:SS(1)
            for j=1:SS(2)                 
                dif_I=(M_I-Im_I(i,j))^2;
                SM_I=SM_I+dif_I;                    
            end
        end
        rms_proposed(counter)=(SM_I/(SS(1)*SS(2)))^0.5; 
        ent_proposed(counter)=entropy(uint16(ImOut)); %Entropy   

    end
end

%Mean rms and entropy
mean_rms_original=mean(rms_original);
mean_rms_proposed=mean(rms_proposed);

mean_entropy_original=mean(ent_original);
mean_entropy_proposed=mean(ent_proposed);

save mean_rms_proposed_Medical mean_rms_proposed
save mean_entropy_proposed_Medical mean_entropy_proposed
save mean_rms_original_Medical mean_rms_original
save mean_entropy_original_Medical mean_entropy_original







          




