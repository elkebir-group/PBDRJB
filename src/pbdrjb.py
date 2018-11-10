import os
import argparse
import subprocess
import numpy as np

class LPSolver(object):
	def __init__(self, args=None, B=None, u=None, e=None, d=None, 
				Z=None, C_ans=None, C_ans_sgn=None, eps=1e-5, alpha=0.99,
				save_dir=None, prob_name=None):
		self.args = args
		self.backend = self.get_lp()
		self.B = B
		self.u = u
		self.e = e
		self.d = d
		self.Z = Z
		self.C_ans = C_ans
		self.C_ans_sgn = C_ans_sgn
		self.eps = eps
		self.alpha = alpha
		self.save_dir = save_dir
		self.prob_name = prob_name

	def get_lp(self):
		from milp import milp
		lp = milp
		return lp

	def get_lp_results(self):
		solution = self.model.solution
		obj_val = solution.get_objective_value()   
		C = np.array([[solution.get_values('c_%s_%s'%(i, p)) \
				for p in range(self.d.shape[0])] for i in range(self.u.shape[0])])
		C_delta = dict()        
		for symbol in ['minus', 'zero', 'plus']:
			C_delta[symbol] = np.array([[solution.get_values('c_%s_%s_%s'%(symbol, i, p)) \
					for p in range(self.d.shape[0])] for i in range(self.u.shape[0])])
			C_delta[symbol] = (C_delta[symbol] > 0.5).astype(int)
		C_sign = (C_delta['plus'] - C_delta['minus']).astype(int)
		return obj_val, C, C_delta, C_sign
	
	def save_solution(self):
		subprocess.call(['mkdir', '-p', self.save_dir])
		fout = open('%s/output_summary.txt'%(self.save_dir), 'w')
		fout.write('Input:\n')
		fout.write('B = \n%s\n'%(self.B))
		fout.write('u = %s\n'%(self.u))
		fout.write('d = %s\n'%(self.d))
		fout.write('e = %s\n'%(self.e))
		if self.C_ans is not None:
			fout.write('C = \n%s\n'%(self.C_ans))	
		if self.C_ans_sgn is not None:
			fout.write('C_sgn = \n%s\n'%(self.C_ans_sgn.astype(int)))
		fout.write('eps = %s\n'%(self.eps))
		fout.write('alpha = %s\n'%(self.alpha))
		fout.write('\nSolution:\n')
		fout.write('Objective value = %s\n'%(self.obj_val))
		fout.write('solution = \n%s\n'%(np.array2string(self.C, suppress_small=True)))
		fout.write('solution_sgn = \n%s\n'%(self.C_sign))
		if self.C_ans is not None:			
			fout.write('norm = %s\n'%(np.linalg.norm(self.C - self.C_ans) / self.C.size))
		fout.close()
		np.savetxt('%s/output_C.txt'%(self.save_dir), self.C, fmt='%g')
		self.model.write('%s/output_model.lp'%(self.save_dir))
	
	def solve(self):
		self.model = self.backend(B=self.B, u=self.u, e=self.e, d=self.d, 
							Z=self.Z, eps=self.eps, alpha=self.alpha)
		self.model.set_results_stream(None)
		self.model.set_warning_stream(None)		
		self.model.solve()
		self.obj_val, self.C, _, self.C_sign = self.get_lp_results()
		self.save_solution()

def parse_input(input_dir):
	B = np.loadtxt('%s/B.txt'%(input_dir))
	u = np.loadtxt('%s/u.txt'%(input_dir))
	e = np.loadtxt('%s/e.txt'%(input_dir))
	d = np.loadtxt('%s/d.txt'%(input_dir))
	if os.path.isfile('%s/C.txt'%(input_dir)):
		C_ans = np.loadtxt('%s/C.txt'%(input_dir))
	else:
		C_ans = None
	if os.path.isfile('%s/C_sign.txt'%(input_dir)):
		C_ans_sgn = np.loadtxt('%s/C_sign.txt'%(input_dir))
	else:
		C_ans_sgn = None
	Z_minus = np.loadtxt('%s/Z_minus.txt'%(input_dir))
	Z_zero = np.loadtxt('%s/Z_zero.txt'%(input_dir))
	Z_plus = np.loadtxt('%s/Z_plus.txt'%(input_dir))
	Z = [Z_minus, Z_zero, Z_plus]
	return B, u, e, d, Z, C_ans, C_ans_sgn

def solve_instance():		
	B, u, e, d, Z, C_ans, C_ans_sgn = parse_input(args.input_dir)	
	solver = LPSolver(args=args, B=B, u=u, e=e, d=d, 
			Z=Z, C_ans=C_ans, C_ans_sgn=C_ans_sgn, eps=args.eps, 
			alpha=args.alpha, save_dir=args.output_dir, prob_name='PDBRJB')
	solver.solve()	

def main():
	solve_instance()
	
if __name__ == '__main__':
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)						
	parser.add_argument('--input_dir', action='store', help='directory of input directory')
	parser.add_argument('--output_dir', action='store', help='directory of output directory')
	parser.add_argument('--alpha', action='store', type=float, default=0.99,
						help='edges with neutral probability greater than this threshold will be removed')
	parser.add_argument('--eps', action='store', type=float, default=1e-5,
						help='neglectable constant used to avoid numerical issues')
	args = parser.parse_args()
	main()
