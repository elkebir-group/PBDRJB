import numpy as np 
import cplex

def milp(B=None, e=None, u=None, Z=None, d=None, eps=1e-5, alpha=0.99):
	"""
	Parameters
	----------
	B : (n, m). Mutation matrix, in which $B_{i,p}$ is 1 if gene $p$ is mutated in clone $i$ and 0 otherwise.
	e : (m, ). Relative abundance vector of the matched normal sample.
	u : (n, ). mixture vector.
	Z : [Z^{-}, Z^{0}, Z^{+}]. Output of xseq. $\mathcal{Z}$ consists of three components $\mathcal{Z}=[Z^{-},Z^{0},Z^{+}]$ 
		where $Z^{-}=[z^{-}_{q,p}], Z^{0}=[z^{0}_{q,p}], Z^{+}=[z^{+}_{q,p}]$. For each sample, 
		the xseq tool gives for each pair of source gene $q$ and target gene $p$ a probability 
		that the target gene $p$ will be down-regulated, neutral or up-regulated due to the mutations 
		in the source gene $q$, which are denoted by $z^{-}_{q,p}, z^{0}_{q,p}, z^{+}_{q,p}$, respectively.
		Z^{-}, Z^{0}, Z^{+} have shape (m, m).
	d : (m, ). Relative abundance vector of the tumor sample.
	eps : a pre-specified constant to convert strict inequalities to inequalities with equality conditions.
    alpha: xseq edges with neutral prob. greater than this threshold will be removed
	
	Returns
	----------
	model : a Cplex model, which will solve the following variables
			C : (n, m). expression matrix.     
	"""
	# n: number of clones; m: number of genes
	n, m = u.shape[0], d.shape[0]

	# unpack Z to Z_minus, Z_zero, Z_plus
	Z_minus, Z_zero, Z_plus = Z
	
	# initialize cplex model
	model = cplex.Cplex()
	
	# maximize the objective function
	model.objective.set_sense(model.objective.sense.maximize)
	
	# set variable: c[i,p]
	model.variables.add(
		lb = [0.0 for _ in range(n) for _ in range(m)],
		ub = [1.0 for _ in range(n) for _ in range(m)],
		names = ['c_%s_%s'%(i, p) for i in range(n) for p in range(m)]
	)
	
	# get log(Z) for Z_minus, Z_zero, Z_plus, this will be the coef. in obj.
	logp = {
		'minus': np.log(Z_minus),
		'zero': np.log(Z_zero),
		'plus': np.log(Z_plus)        
	}              
	
	E = (Z_zero < alpha - eps)
	in_q = [[np.where(B[i, :] * E[:, p] == 1)[0] for p in range(m)] for i in range(n)]
	
	for symbol in ['minus', 'zero', 'plus']:
		coef = np.log(alpha) if symbol == 'zero' else np.log((1 - alpha) / 2.)
		model.variables.add(
			obj = [1.0 * (np.sum(logp[symbol][in_q[i][p], p]) + (len(in_q[i][p]) == 0) * coef)
					for i in range(n) for p in range(m)],
			lb = [0 for _ in range(n) for _ in range(m)],
			ub = [1 for _ in range(n) for _ in range(m)],
			types = ['B' for i in range(n) for p in range(m)],
			names = ['c_%s_%s_%s'%(symbol, i, p) for i in range(n) for p in range(m)]
		)
	# add constraints: every row sum is 1 in matrix C    
	model.linear_constraints.add(
		lin_expr = [[['c_%s_%s'%(i, p) for p in range(m)], np.ones(m)] for i in range(n)],
		rhs = [1.0 - eps] * n,
		senses = ['G'] * n
	)
	model.linear_constraints.add(
		lin_expr = [[['c_%s_%s'%(i, p) for p in range(m)], np.ones(m)] for i in range(n)],
		rhs = [1.0 + eps] * n,
		senses = ['L'] * n
	)
	   
	# add constraints: d=uC
	model.linear_constraints.add(
		lin_expr = [[['c_%s_%s'%(i, p) for i in range(n)], u] for p in range(m)],
		rhs = d,
		senses = ['E'] * m
	)
	
	# add constraints: exatly one of c_minus[i,p], c_zero[i,p], c_plus[i,p] is 1
	model.linear_constraints.add(
		lin_expr = [[['c_%s_%s_%s'%(symbol, i, p) for symbol in ['minus', 'zero', 'plus']], 
					 [1, 1, 1]]
				   for i in range(n) for p in range(m)],
		rhs = [1] * n * m,
		senses = ['E'] * n * m
	)

	# add constraints: c[i,p] is down-regulated (c_minus[i,p] is 1) 
	model.linear_constraints.add(
		lin_expr = [[['c_minus_%s_%s'%(i, p), 'c_%s_%s'%(i, p)], [1, 1]]
				   for i in range(n) for p in range(m)],
		rhs = [e[p] for i in range(n) for p in range(m)],
		senses = ['G'] * n * m
	)


	# add constraints: c[i,p] is up-regulated (c_plus[i,p] is 1)
	model.linear_constraints.add(
		lin_expr = [[['c_plus_%s_%s'%(i, p), 'c_%s_%s'%(i, p)], [1, -1]]
				   for i in range(n) for p in range(m)],
		rhs = [-e[p] for i in range(n) for p in range(m)],
		senses = ['G'] * n * m
	)
	
	# add aux variable to model multiplication
	# f_minus[i,p] = c_minus[i,p] * c[i,p]
	# f_plus[i,p] = c_plus[i,p] * c[i,p]
	for symbol in ['minus', 'plus']:
		model.variables.add(            
			lb = [0 for _ in range(n) for _ in range(m)],
			ub = [1 for _ in range(n) for _ in range(m)],            
			names = ['f_%s_%s_%s'%(symbol, i, p) for i in range(n) for p in range(m)]
		)
		# add constraints to model multiplication {0, 1} x [0, 1]: 
		model.linear_constraints.add(
			lin_expr = [[['f_%s_%s_%s'%(symbol, i, p), 'c_%s_%s_%s'%(symbol, i, p)], [1, -1]] 
						for i in range(n) for p in range(m)],
			rhs = [0] * n * m,
			senses = ['L'] * n * m
		)
		model.linear_constraints.add(
			lin_expr = [[['f_%s_%s_%s'%(symbol, i, p), 'c_%s_%s'%(i, p)], [1, -1]] 
						for i in range(n) for p in range(m)],
			rhs = [0] * n * m,
			senses = ['L'] * n * m
		)
		model.linear_constraints.add(
			lin_expr = [[['f_%s_%s_%s'%(symbol, i, p), 'c_%s_%s_%s'%(symbol, i, p), 'c_%s_%s'%(i, p)], [1, -1, -1]] 
						for i in range(n) for p in range(m)],
			rhs = [-1] * n * m,
			senses = ['G'] * n * m
		)

	# add constraints: c_minus[i,p] * c[i,p]  - (1 - c_minus[i,p]) <= max{e[p] - eps, 0}
	model.linear_constraints.add(
		lin_expr = [[['f_minus_%s_%s'%(i,p), 'c_minus_%s_%s'%(i, p)], [1, 1]] for i in range(n) for p in range(m)],
		rhs = [max(e[p] - eps, 0) + 1 for i in range(n) for p in range(m)],
		senses = ['L'] * n * m
	)
	model.linear_constraints.add(
		lin_expr = [[['f_plus_%s_%s'%(i,p), 'c_plus_%s_%s'%(i,p)], [1, -1]] for i in range(n) for p in range(m)],
		rhs = [min(e[p] + eps, 1) - 1 for i in range(n) for p in range(m)],
		senses = ['G'] * n * m
	)

	return model