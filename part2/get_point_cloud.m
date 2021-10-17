function [pc] = get_point_cloud(rgb_image, depth_array, cam_params, scale_factor)
    
    %%Generate XYZ from depth
    I = (1:480)'*ones(1,640);
    J = ones(480,1)*(1:640);
    one = ones(480,640);

    pixel_map = [J(:)'; I(:)'; one(:)'];
    d = double(depth_array(:)')/scale_factor; %convert to meters
    triple_depth = [d; d; d];

    xyz = (cam_params.Kdepth)\(triple_depth.*pixel_map);

    %Rigid transformation - RGB camera
    xyz_rgb = cam_params.R  * xyz + cam_params.T*ones(1, size(xyz,2));
    %Project in the RGB image plane (homog. coordinates)
    omega = cam_params.Krgb * xyz_rgb;

    %Euclidean coordinates
    u = round(omega(1.,:)./omega(3,:));
    v = round(omega(2,:)./omega(3,:));

    %just clean up wrong coordinates
    v(v > 480) = 1;
    v(v < 1) = 1;
    u(u > 640) = 1;
    u(u < 1) = 1; %%
    
    %convert to linear index
    ind = sub2ind(size(depth_array),v,u);

    %Generate "virtual" depth image in RGB
    depth_new = zeros(size(depth_array));
    depth_new(ind) = xyz_rgb(3,:);

    %Get fully "colored" point cloud
    imnew = reshape(rgb_image,[640* 480 3]);
    imnew(depth_new(:) == 0,:) = 0;
    
    xyz_rgb(1,sum(imnew,2) == 0) = 0;
    xyz_rgb(2,sum(imnew,2) == 0) = 0;
    xyz_rgb(3,sum(imnew,2) == 0) = 0;
    pc = pointCloud(xyz_rgb','color',imnew);

end
