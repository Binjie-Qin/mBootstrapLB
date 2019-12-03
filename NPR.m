function B = NPR(X,Y,h)
sizeX = size(X);
sizeY = size(Y);
A = zeros(3,3);
b = zeros(3,1);
beta = zeros(3,sizeX(1));
temp_reg = zeros(sizeX(1),1);
%k = 0;
for i = 3:(sizeX(1)-2)
    if i <= (h/2)
        L = i - 1;
        for j = -L:L
            kernel_weight = 1-log((exp(1)-1)*abs(j/L)+1);
            b(1) = b(1) + Y(i+j)*kernel_weight;
            b(2) = b(2) + Y(i+j)*j*kernel_weight;
            b(3) = b(3) + Y(i+j)*j^2*kernel_weight;
            A(1,1) = A(1,1)+1*kernel_weight;
            A(1,2) = A(1,2)+j*kernel_weight;
            A(1,3) = A(1,3)+j^2*kernel_weight;
            A(2,1) = A(2,1)+j*kernel_weight;
            A(2,2) = A(2,2)+j^2*kernel_weight;
            A(2,3) = A(2,3)+j^3*kernel_weight;
            A(3,1) = A(3,1)+j^2*kernel_weight;
            A(3,2) = A(3,2)+j^3*kernel_weight;
            A(3,3) = A(3,3)+j^4*kernel_weight;
        end
        %beta(:,i) = A^(-1)*b;
        %temp_reg(i) = beta(1,i);
        if rank(A) == 3
            beta(:,i) = A^(-1)*b;
            temp_reg(i) = beta(1,i);
        else
            temp_reg(i) = Y(i);
        end
        A = zeros(3,3);
        b = zeros(3,1);
    else if i >= (sizeX(1)+1-h/2)
            L = sizeX(1) - i;
            for j = -L:L
                kernel_weight = 1-log((exp(1)-1)*abs(j/L)+1);
                b(1) = b(1) + Y(i+j)*kernel_weight;
                b(2) = b(2) + Y(i+j)*j*kernel_weight;
                b(3) = b(3) + Y(i+j)*j^2*kernel_weight;
                A(1,1) = A(1,1)+1*kernel_weight;
                A(1,2) = A(1,2)+j*kernel_weight;
                A(1,3) = A(1,3)+j^2*kernel_weight;
                A(2,1) = A(2,1)+j*kernel_weight;
                A(2,2) = A(2,2)+j^2*kernel_weight;
                A(2,3) = A(2,3)+j^3*kernel_weight;
                A(3,1) = A(3,1)+j^2*kernel_weight;
                A(3,2) = A(3,2)+j^3*kernel_weight;
                A(3,3) = A(3,3)+j^4*kernel_weight;
            end
            %beta(:,i) = A^(-1)*b;
            %temp_reg(i) = beta(1,i);
            if rank(A) == 3
                beta(:,i) = A^(-1)*b;
                temp_reg(i) = beta(1,i);
            else
                temp_reg(i) = Y(i);
            end
            A = zeros(3,3);
            b = zeros(3,1);
        else
            L = h/2;
            for j = -L:L
                kernel_weight = 1-log((exp(1)-1)*abs(j/L)+1);
                b(1) = b(1) + Y(i+j)*kernel_weight;
                b(2) = b(2) + Y(i+j)*j*kernel_weight;
                b(3) = b(3) + Y(i+j)*j^2*kernel_weight;
                A(1,1) = A(1,1)+1*kernel_weight;
                A(1,2) = A(1,2)+j*kernel_weight;
                A(1,3) = A(1,3)+j^2*kernel_weight;
                A(2,1) = A(2,1)+j*kernel_weight;
                A(2,2) = A(2,2)+j^2*kernel_weight;
                A(2,3) = A(2,3)+j^3*kernel_weight;
                A(3,1) = A(3,1)+j^2*kernel_weight;
                A(3,2) = A(3,2)+j^3*kernel_weight;
                A(3,3) = A(3,3)+j^4*kernel_weight;
            end
            beta(:,i) = A^(-1)*b;
            temp_reg(i) = beta(1,i);
            if rank(A) == 3
                beta(:,i) = A^(-1)*b;
                temp_reg(i) = beta(1,i);
            else
                temp_reg(i) = Y(i);
            end
            A = zeros(3,3);
            b = zeros(3,1);
        end
    end
end
temp_reg(1:2) = Y(1:2);
temp_reg(sizeX(1)-1:sizeX(1)) = Y(sizeX(1)-1:sizeX(1));
B = temp_reg;
k = 0;
while k < 10
    k = k + 1;
    for i = 3:(sizeX(1)-2)
        if i <= (h/2)
            L = i - 1;
            res = abs(Y((i-L):(i+L)) - B((i-L):(i+L)));
            sigma = median(res)/0.6745+eps;
            for j = -L:L
                robust_weight = max((1-((Y(i+j)-B(i+j))/(3*sigma))^2),0);
                robust_weight = robust_weight^2;
                kernel_weight = 1-log((exp(1)-1)*abs(j/L)+1);
                b(1) = b(1) + Y(i+j)*kernel_weight*robust_weight;
                b(2) = b(2) + Y(i+j)*j*kernel_weight*robust_weight;
                b(3) = b(3) + Y(i+j)*j^2*kernel_weight*robust_weight;
                A(1,1) = A(1,1)+1*kernel_weight*robust_weight;
                A(1,2) = A(1,2)+j*kernel_weight*robust_weight;
                A(1,3) = A(1,3)+j^2*kernel_weight*robust_weight;
                A(2,1) = A(2,1)+j*kernel_weight*robust_weight;
                A(2,2) = A(2,2)+j^2*kernel_weight*robust_weight;
                A(2,3) = A(2,3)+j^3*kernel_weight*robust_weight;
                A(3,1) = A(3,1)+j^2*kernel_weight*robust_weight;
                A(3,2) = A(3,2)+j^3*kernel_weight*robust_weight;
                A(3,3) = A(3,3)+j^4*kernel_weight*robust_weight;
            end
            if rank(A) == 3
                beta(:,i) = A^(-1)*b;
                temp_reg(i) = beta(1,i);
            else
                temp_reg(i) = B(i);
            end
            A = zeros(3,3);
            b = zeros(3,1);
            clear res;
        else if i >= (sizeX(1)+1-h/2)
                L = sizeX(1) - i;
                res = abs(Y((i-L):(i+L)) - B((i-L):(i+L)));
                sigma = median(res)/0.6745+eps;
                for j = -L:L
                    robust_weight = max((1-((Y(i+j)-B(i+j))/(3*sigma))^2),0);
                    robust_weight = robust_weight^2;
                    kernel_weight = 1-log((exp(1)-1)*abs(j/L)+1);
                    b(1) = b(1) + Y(i+j)*kernel_weight*robust_weight;
                    b(2) = b(2) + Y(i+j)*j*kernel_weight*robust_weight;
                    b(3) = b(3) + Y(i+j)*j^2*kernel_weight*robust_weight;
                    A(1,1) = A(1,1)+1*kernel_weight*robust_weight;
                    A(1,2) = A(1,2)+j*kernel_weight*robust_weight;
                    A(1,3) = A(1,3)+j^2*kernel_weight*robust_weight;
                    A(2,1) = A(2,1)+j*kernel_weight*robust_weight;
                    A(2,2) = A(2,2)+j^2*kernel_weight*robust_weight;
                    A(2,3) = A(2,3)+j^3*kernel_weight*robust_weight;
                    A(3,1) = A(3,1)+j^2*kernel_weight*robust_weight;
                    A(3,2) = A(3,2)+j^3*kernel_weight*robust_weight;
                    A(3,3) = A(3,3)+j^4*kernel_weight*robust_weight;
                end
                %beta(:,i) = A^(-1)*b;
                %temp_reg(i) = beta(1,i);
                if rank(A) == 3
                    beta(:,i) = A^(-1)*b;
                    temp_reg(i) = beta(1,i);
                else
                    temp_reg(i) = B(i);
                end
                A = zeros(3,3);
                b = zeros(3,1);
                clear res;
            else
                L = h/2;
                res = abs(Y((i-L):(i+L)) - B((i-L):(i+L)));
                sigma = median(res)/0.6745+eps;
                for j = -L:L
                    robust_weight = max((1-((Y(i+j)-B(i+j))/(3*sigma))^2),0);
                    robust_weight = robust_weight^2;
                    kernel_weight = 1-log((exp(1)-1)*abs(j/L)+1);
                    b(1) = b(1) + Y(i+j)*kernel_weight*robust_weight;
                    b(2) = b(2) + Y(i+j)*j*kernel_weight*robust_weight;
                    b(3) = b(3) + Y(i+j)*j^2*kernel_weight*robust_weight;
                    A(1,1) = A(1,1)+1*kernel_weight*robust_weight;
                    A(1,2) = A(1,2)+j*kernel_weight*robust_weight;
                    A(1,3) = A(1,3)+j^2*kernel_weight*robust_weight;
                    A(2,1) = A(2,1)+j*kernel_weight*robust_weight;
                    A(2,2) = A(2,2)+j^2*kernel_weight*robust_weight;
                    A(2,3) = A(2,3)+j^3*kernel_weight*robust_weight;
                    A(3,1) = A(3,1)+j^2*kernel_weight*robust_weight;
                    A(3,2) = A(3,2)+j^3*kernel_weight*robust_weight;
                    A(3,3) = A(3,3)+j^4*kernel_weight*robust_weight;
                end
                %beta(:,i) = A^(-1)*b;
                %temp_reg(i) = beta(1,i);
                if rank(A) == 3
                    beta(:,i) = A^(-1)*b;
                    temp_reg(i) = beta(1,i);
                else
                    temp_reg(i) = B(i);
                end
                A = zeros(3,3);
                b = zeros(3,1);
                clear res;
            end
        end
    end
    k = k
    %temp_reg(1:2) = Y(1:2);
    temp_reg(1:2) = 0;
    %temp_reg(2) = beta(1,3) - beta(2,3) + beta(3,3);
    %temp_reg(1) = beta(1,3) - 2*beta(2,3) + 4*beta(3,3);
    temp_reg(sizeX(1)-1:sizeX(1)) = 0;
    %temp_reg(sizeX(1)-1:sizeX(1)) = Y(sizeX(1)-1:sizeX(1));
    %temp_reg(sizeX(1)-1) = beta(1,sizeX(1)-2) + beta(2,sizeX(1)-2) + beta(3,sizeX(1)-2);
    %temp_reg(sizeX(1)) = beta(1,sizeX(1)-2) + 2*beta(2,sizeX(1)-2) + 4*beta(3,sizeX(1)-2);
    B = temp_reg;
end

%     temp_reg(2) = beta(1,3) - beta(2,3) + beta(3,3);
%     temp_reg(1) = beta(1,3) - 2*beta(2,3) + 4*beta(3,3);
%     temp_reg(sizeX(1)-1) = beta(1,sizeX(1)-2) + beta(2,sizeX(1)-2) + beta(3,sizeX(1)-2);
%     temp_reg(sizeX(1)) = beta(1,sizeX(1)-2) + 2*beta(2,sizeX(1)-2) + 4*beta(3,sizeX(1)-2);
    B = temp_reg;
end