def PrintBanner(eig_solver, verbosity=1):
    if verbosity>=1:
        if verbosity==2:
            print('\nSDPSubspaceSolver called with eigenvalue solver', eig_solver, 'and SDP solver cvx')
            print('    |                                         | Runtime Ratios | ')
            print(' it | probj        dures    prres    dugap    | eig  SDP  main | eigtol   prdim  time |  y[0]         y[1]')
        elif verbosity==1:
            print('\nSDPSubspaceSolver called with eigenvalue solver', eig_solver, 'and SDP solver cvx')
            print('    |                                         | Runtime Ratios | ')
            print(' it | probj        dures    prres    dugap    | eig  SDP  main | eigtol   prdim  time')
    return