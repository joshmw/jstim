function makegabors

for i = 1:100
grating = mglMakeGrating(10,10,1.8,round(rand*360),0);
name = ['image' num2str(i) '.png']
imwrite(grating,name)
end




