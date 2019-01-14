i=2;
j=2;
val=0;
while(j<= size(image,2)-1);
    while(i<= size(image,1)-1); 
       val=image(i,j);
       if(val ~= 0);
           if (image(i+1,j)<val*10);
               image(i+1,j)=val*10;   %% don't forget the side cases which might have problems
           end;
           if (image(i,j+1)<val*10);
               image(i,j+1)=val*10;
           end;
           if (image(i-1,j)<val*10);
               image(i-1,j)=val*10;
           end;
           if(image(i,j-1)<val*10);
           image(i,j-1)=val*10; %% i might want to have the corners also included
           end;
           image(i,j)=val*15;
       end;
       i=i+1;
    end;
    j=j+1;
end;