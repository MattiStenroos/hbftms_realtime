function [coil] = make_coil_moment(index_xy, index_z)
%MAKE_COIL_MOMENT makes a coil model using precomputed moments. 
%
%   [coil] = MAKE_COIL_MOMENT(index_xy, index_z)
%
%   index_xy = 1 ... 9
%       correspond to the coils mentioned in (Stenroos and Koponen, 2019).
%   index xy = 10 ... 13
%       correspond to next four moment coils found with numerical methods.
%
%   index_z = 1 ... 3
%       number of layers
%
%   The coils have PPL dipoles per layer, where PPL follows from index_xy:
%   index_xy   PPL
%          1     2
%          2     6
%          3     8
%          4    12
%          5    14
%          6    22
%          7    26
%          8    38
%          9    42
%         10    56
%         11    62
%         12    80
%         13    98
%
%   v191002 (c) Lari Koponen (lari.koponen@aalto.fi

    % Load parameters for moment coils.
    [xy, z] = parameters();
       
    % Stack thin moment coils as described by z moments
    coil = [];
    coil.QP = [];
    coil.QN = []; 
    coil.QW = [];
    [pts, wts] = generate_height(z(index_z, :));
    for i = 1 : length(pts)
        layer = generate_layer(xy(index_xy, :), i);
        coil.QP = [coil.QP; bsxfun(@plus,layer.QP,[0 0 1] * pts(i))];
        coil.QN = [coil.QN; layer.QN];
        coil.QW = [coil.QW; wts(i) * layer.QW];
    end
    coil.QPinds = [1 size(coil.QP, 1)];
    
end

function [pts, wts] = generate_height(row)
%GENERATE_HEIGHT Generates placement of layers of a figure-of-eight coil.
%
%   [pts, wts] = GENERATE_HEIGHT(row)
%
% Version 2019-09-12

    pts = [];
    wts = [];
    if mod(row(2), 2) == 1
        pts = [pts; 0];
        wts = [wts; row(4)];
    end
    if ~isnan(row(3))
        pts = [pts; -row(3); row(3)];
        wts = [wts; row(5); row(5)];
    end
    
    % Scale to correct thickness (of windings and plastic casing)
    pts = (0.003 + 0.5 * 0.007) + 0.007 * pts;
     
    % Sort from bottom up
    [~, permutation] = sort(pts);
    pts = pts(permutation);
    wts = wts(permutation);
end

function [coil] = generate_layer(row, index)
%GENERATE_LAYER Generates one layer of figure-of-eight coil.
%
%   [coil] = GENERATE_LAYER(row, index)
%
% where row is a row of input data and index is the z-index of the layer.
%
% Version 2019-10-02

    pts = [];
    wts = [];
    
    % For all coils, the adjacent rings must have half-a-dipole offset.
    % And, for all coils, the outermost rings must have their dipoles as
    % far from the origin as possible. The two are controlled by "offset".
    if row(5) > 0 || row(4) == 0
        offset = mod(index, 2);
    else
        offset = mod(index + 1, 2);
    end

    % Add point to origin (if any)
    if row(2)
        pts = [pts; 0 0];
        wts = [wts; row(9)];
    end
    % Create innermost ring (if any)
    if row(3)
        theta = linspace(0, 2 * pi, row(3) + 1);
        if offset
            theta = theta + 0.5 * (theta(2) - theta(1));
        end
        theta = theta(1 : end - 1).';
        temp = row(6) * [cos(theta) sin(theta)];
        pts = [pts; temp];
        wts = [wts; ones(row(3), 1) * row(10)];
    end
    % Create second ring (if any), offset it from first ring
    if row(4)
        theta = linspace(0, 2 * pi, row(4) + 1);
        if ~offset
            theta = theta + 0.5 * (theta(2) - theta(1));
        end
        theta = theta(1 : end - 1).';
        temp = row(7) * [cos(theta) sin(theta)];
        pts = [pts; temp];
        wts = [wts; ones(row(4), 1) * row(11)];
    end
    % Create third ring (if any), offset it from second ring
    if row(5)
        theta = linspace(0, 2 * pi, row(5) + 1);
        if offset
            theta = theta + 0.5 * (theta(2) - theta(1));
        end
        theta = theta(1 : end - 1).';
        temp = row(8) * [cos(theta) sin(theta)];
        pts = [pts; temp];
        wts = [wts; ones(row(5), 1) * row(12)];
    end
    
    % Build figure-of-eight coil from two such circular coils
    pts(:, 1) = 0.044 - pts(:, 1);
    pts(:, 3) = 0;
    
    pts = [pts; bsxfun(@times,[-1 1 1],pts)];
    wts = [wts; -wts];
    
    coil = [];
    coil.QP = pts;
    coil.QN = ones(size(pts, 1), 1) * [0 0 1];
    coil.QW = wts;

end

function [xy, z] = parameters()
%PARAMETERS Returns the parameters for moment coils.
%   More specifically, the function returns xy moments for
%   Magstim 70mm Double Coil, based on the dimensions given by
%   (Thielscher and Kammer, 2002). The returned parameters for z are 
%   universal to any coil.
%
%   For computation of the moment-coil parameters mentioned in
%   (Stenroos and Koponen, 2019), see attached Mathematica notebook.
%
%   [xy, z] = PARAMETERS()
%
% Version 2019-09-12

    xy = [ ...
        1,  1  0  0  0, nan                     nan                     nan,                     0.035489566897954745481  nan                       nan                       nan
        2,  0  3  0  0, 0.026206805412127250834 nan                     nan,                     nan                      0.011829855632651581827   nan                       nan
        3,  0  4  0  0, 0.026206805412127250834 nan                     nan,                     nan                      0.0088723917244886863703  nan                       nan
        4,  1  5  0  0, 0.031386054575977978885 nan                     nan,                     0.010746392494122686594  0.0049486348807664117775  nan                       nan
        5,  1  6  0  0, 0.031386054575977978885 nan                     nan,                     0.010746392494122686594  0.0041238624006386764812  nan                       nan
        6,  1  5  5  0, 0.027582967775204804917 0.042470242610439607606 nan,                     0.0074361826276453370784 0.0050295026949271764727  0.00058117415913470520778 nan
        7,  1  6  6  0, 0.028572123052414024716 0.045532807828329811726 nan,                     0.0080624506066454296814 0.0043081591207255534349  0.00026302692782599919834 nan
        8,  1  9  9  0, 0.023805440173486153010 0.037734505536486966404 nan,                     0.0051167963914743284062 0.0024464134652038830332  0.00092833881329394108627 nan
        9,  1 10 10  0, 0.023805440173486153010 0.037734505536486966404 nan,                     0.0051167963914743284062 0.0022017721186834947298  0.00083550493196454697764 nan
        10, 1  9  9  9, 0.020814419673387527712 0.033752508145072085053 0.042002403207437547747, 0.0037239292766285363854 0.0020203332946784583161  0.0013273531126072773712  0.00018182888397273198995
        11, 1 10 10 10, 0.021069102417414577213 0.034173460971248848754 0.042500825022381031525, 0.0038216049680111819857 0.0018575839519868048425  0.0011778320591186816379  0.00013138018188886986918
        12, 1 13 13 13, 0.018959174664545973634 0.031346764571316691615 0.040351016925443498955, 0.0030078506259804427947 0.0012007390031263454506  0.0010375811509768556075  0.00026027340527943761015
        13, 1 16 16 16, 0.018959174664545973634 0.031346764571316691615 0.040351016925443498955, 0.0030078506259804427947 0.00097560044004015567858 0.00084303468516869518108 0.00021147214178954305824
        ];

    z = [ ...
        1, 1, nan,          1   nan
        3, 2, 1 / sqrt(12), nan 1/2
        5, 3, sqrt(3 / 20), 4/9 5/18
        ];
    
end