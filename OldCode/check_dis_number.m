function count=check_dis_number(y)

y = y(:);

count = 1; 
list = y(1);

for i = 2:length(y)
    
    if sum(y(i) == list)==0

        count = count + 1;
        list = [list; y(i)];
        
    end
       
    
end