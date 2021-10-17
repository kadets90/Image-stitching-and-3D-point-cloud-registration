clear
load sintetic1
%load sinteticotable
%load sinteticoprojective
%load sinteticoaffine
%load sinteticoaffine2

%image_list = image_list(1:9)';
[ Homog, world_pts_x, world_pts_y]=part1(image_list, match_list );

minx=min(world_pts_x(:));miny=min(world_pts_y(:));
maxx=max(world_pts_x(:));maxy=max(world_pts_y(:));
mosaico=zeros(fix(maxy-miny)+1,fix(maxx-minx)+1,3);
mosaico2=mosaico;
axis ij;axis equal;
%just for fun
for i=1:length(image_list),
    im2=imread(image_list{i});
    tt=projective2d(Homog{i,1}');
    im3c1=imwarp(uint8(ones(size(im2))),tt,'OutputView',imref2d([size(mosaico,1) size(mosaico,2)],[minx maxx],[miny maxy]));
    im3c=imwarp(im2,tt,'OutputView',imref2d([size(mosaico,1) size(mosaico,2)],[minx maxx],[miny maxy]));
    mosaico=mosaico+double(im3c);
    mosaico2=mosaico2+double(im3c1);
    figure(1);plot([world_pts_x(i,:) world_pts_x(i,1)]',[world_pts_y(i,:) world_pts_y(i,1)]');   
    axis ij;axis equal; hold on;
    figure(2);imagesc(im3c);
    figure(3);imagesc((mosaico./mosaico2)/255);
    drawnow;
end
