function [xyz_rgb] = get_3D_coordinates(idx, depth_array, cam_params, scale_factor)
        
    % Pass 2D coordinates to 3D
    %%Generate XYZ from depth
    one = ones(1,size(idx,1));

    pixel_map = [idx'; one];
    indices = sub2ind(size(depth_array), idx(:,2)', idx(:,1)');
    depth_array = depth_array(indices);
    d = double(depth_array(:)')/scale_factor; %convert to meters
    triple_depth = [d; d; d];

    xyz = (cam_params.Kdepth)\(triple_depth.*pixel_map);

    %Rigid transformation - RGB camera
    xyz_rgb = cam_params.R  * xyz + cam_params.T*ones(1, size(xyz,2));
    xyz_rgb = xyz_rgb';
end