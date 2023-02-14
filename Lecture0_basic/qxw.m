% qxw.m

COUNT = 10000000;
count = 0;

answer = [1,2,3,4,5]; % 答案

for i = 1:COUNT
    
    youranswer = randperm(7);
    youranswer = youranswer(1:5); % 你蒙的答案
    
    right = 0; % 有没有对的
    for j = 1:5
        if youranswer(j) == answer(j) % 对了一个
            right = 1;
        end
    end
    
    if right == 0 % 没有对的
        count = count + 1;
    end
end

p = count/COUNT;