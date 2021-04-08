def calculate_ripley(radii, sample_size, d1=None, d2=None, d3=None, sample_shape='circle', boundary_correct=False, CSR_Normalise=False):
    results = []
    tree, dimensions = make_tree(d1=d1, d2=d2, d3=d3)
    if type(radii) is not list:
        radii = [radii]
    for radius in radii:
        if dimensions == 1:
            score_vol = radius*2
            bound_size = sample_size
            counts = 0
            for x in zip(d1):
                if boundary_correct:
                    vol = calculate_overlap([x], sample_size, radius, sample_shape, dimensions)
                    boundary_correction = vol/score_vol
                    counts += (len(tree.query_ball_point([x], radius))-1)/boundary_correction
                else:
                    counts += len(tree.query_ball_point([x], radius))-1
        elif dimensions == 2:
            score_vol = np.pi * radius**2
            if sample_shape=='circle':
                bound_size = np.pi * sample_size**2
            elif sample_shape=='rectangle':
                bound_size = sample_size[0]*sample_size[1]
            counts = 0
            for x, y in zip(d1, d2):
                if boundary_correct:
                    vol = calculate_overlap([x,y], sample_size, radius, sample_shape, dimensions)
                    boundary_correction = vol/score_vol
                    counts += (len(tree.query_ball_point([x,y], radius))-1)/boundary_correction
                else:
                    counts += len(tree.query_ball_point([x,y], radius))-1
        else:
            score_vol = (4/3) * np.pi * radius**3
            if sample_shape=='circle':
                bound_size = (4/3) * np.pi * sample_size**3
            elif sample_shape=='rectangle':
                bound_size = sample_size[0]*sample_size[1]*sample_size[2]
            counts = 0
            for x, y, z in zip(d1, d2, d3):
                if boundary_correct:
                    vol = calculate_overlap([x,y,z], sample_size, radius, sample_shape, dimensions)
                    boundary_correction = vol/score_vol
                    counts += (len(tree.query_ball_point([x,y,z], radius))-1)/boundary_correction
                else:
                    counts += len(tree.query_ball_point([x,y,z], radius))-1
        if CSR_Normalise:
            results.append((bound_size*counts/len(d1)**2) - score_vol)
        else:
            results.append(bound_size*counts/len(d1)**2)
    if len(results)==1:
        return results[0]
    else:
        return results
