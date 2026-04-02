clear;clc;close all;
folder = "sounds/";
d_ = dir(folder);
names = string({d_.name}).';
names = names(3:end-1);

names_static = names(contains(names,"static"));
names_static = names_static(~contains(names_static,"(0)"));

soi_int_pairs = [];

for name = names_static.'
    ch_name = char(name);
    names_no_static = names(~contains(names,"static"));
    names_no_static = names_no_static(~contains(names_no_static, ch_name(1)));
    names_no_static = names_no_static(~contains(names_no_static, "x"));

    soi_int_pairs = [soi_int_pairs;[repmat(folder+name,length(names_no_static),1), folder+names_no_static]];
end
save('pairs_aux.mat','soi_int_pairs');
