% qxw.m

COUNT = 10000000;
count = 0;

answer = [1,2,3,4,5]; % ��

for i = 1:COUNT
    
    youranswer = randperm(7);
    youranswer = youranswer(1:5); % ���ɵĴ�
    
    right = 0; % ��û�жԵ�
    for j = 1:5
        if youranswer(j) == answer(j) % ����һ��
            right = 1;
        end
    end
    
    if right == 0 % û�жԵ�
        count = count + 1;
    end
end

p = count/COUNT;