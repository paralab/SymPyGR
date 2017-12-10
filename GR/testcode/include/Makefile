all:  bssn_rhs_files

bssn_rhs_files: gen_code
	./gen_code
	./gen_opt_code
	./gen_vec_code

clean:
	/bin/rm -f New*.h MR_deriv*h opt_*.h vec_*.h
