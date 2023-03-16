import Alignment

def read_Fasta(txt_name, path =''):
	"""
	讀取txt檔中的sequences存成dict
	"""

	#把指針移到結尾(判斷結尾)
	f = open(path + txt_name + '.txt', 'r')	                    
	f.seek(0,2)            
	size=f.tell()           
	f.seek(0,0) 

	speices_dic = {}	
	i = 0

	#用指針判斷是否在結尾前
	while f.tell() < size:
		for line in f.readlines():	
			if len(line) <= 0:
				continue		
			if '>' in line:	
				i += 1
				if i != 1:
					speices_dic[speices_number] = res

				speices_number = line.replace('>','').strip('\n')
				res = ''
			else:
				res = res + str(line).strip()

	speices_dic[speices_number] = res  
	f.close()

	return speices_dic

def isdigit(s):
	if type(s) == int or type(s) == float:
		return True
	elif type(s) == list:
		for i in range(len(s)):
			if type(s[i]) != int and type(s[i]) != float:
				return False
		return True
	else:
		return False

def is_AA(s):
	Amino_acid = ['A','B','C','D','E','F','G','H']
	for a in s:
		if a not in Amino_acid:
			return False
	return True


def get_seqs():
	seqx = input('Please enter seq1(life X) TXT FILE in Fasta form: ')
	try:
		lifex = read_Fasta(seqx, path ='')
	except FileNotFoundError:
		print(f'----------File not found----------\n')
		return get_seqs()
	except:
		print('Need to be FASTA\n')
		return get_seqs()
	seq = lifex[list(lifex.keys())[0]].strip() 
	if not is_AA(seq):
		print(f'----------Amino acid in Life X out of range----------\n')
		return get_seqs()

	ref = input('Please enter reference sequence TXT FILE in Fasta(Could be multyple seqs): ')
	try:
		ref = read_Fasta(ref, path ='')
	except FileNotFoundError:
		print(f'----------File not found----------\n')
		return get_seqs()
	except:
		print(f'----------Need to be FASTA----------\n')
		return get_seqs()

	return 	lifex , ref

def go_align(lifex , ref ):
	try:		
		aligntype = input("Global alignment enter 0, local alignment enter 1, Semi-Global enter 2: ")
		if aligntype not in ['0','1','2']:
			print(f'----------Please enter 0 or 1----------\n')
			return go_align(lifex , ref )
		p = input('If you want to print the result enter (y), or enter (n): ')
		print('----------Enter the paramaters you want in this alginment----------')
		try:
			match_score = float(input(f'Please enter match score: '))
		except ValueError:
			print(f'----------match score most be digit----------\n')	
			return go_align(lifex , ref )
		try:
			mismatch_score = float(input(f'Please enter mismatch score: '))
		except ValueError:
			print(f'----------mismatch score most be digit----------\n')	
			return go_align(lifex , ref )
		try:
			gap_penalty = float(input(f'Please enter gap_penalty: '))
		except ValueError:
			print(f'----------gap_penalty most be digit----------\n')	
			return go_align(lifex , ref )

		affine = input('Affine penalty enter 0, linear penalty enter 1: ')
		if affine != '0' and affine != '1':
			print(f'----------Plase enter 0 or 1----------\n')
			return go_align(lifex , ref )
		elif affine == '0':
			affine = True
			try:
				extend_gap_penalty = float(input('Please enter extend gap penalty: '))
			except ValueError:
				print(f'----------extend gap penalty most be digit----------\n')
				return go_align(lifex , ref )	
		else:
			affine = False
			extend_gap_penalty = gap_penalty
		print('--------------------Alignment Start--------------------\n')
		i = 1
		seq = lifex[list(lifex.keys())[0]].strip() #lifex
		for speices in list(ref.keys()):
			seq_query = ref[speices].strip()
			print(f'<第{i}次Alignemnt...>')
			print('-----------------------')
			if not is_AA(seq_query):
				print(f'----------Amino acid in {i} species out of range----------\n')
				i += 1
				continue
			if aligntype == '0':
				alignment = Alignment.Global_Alignment(seq, seq_query, gap_penalty,extend_gap_penalty, match_score, mismatch_score)
			elif aligntype == '1':
				alignment = Alignment.Local_Alignment(seq, seq_query, gap_penalty,extend_gap_penalty, match_score, mismatch_score)
			elif aligntype == '2':
				alignment = Alignment.SemiGlobal_Alignment(seq, seq_query, gap_penalty,extend_gap_penalty, match_score, mismatch_score)
			alignment.align(affine = affine)
			if p == 'y':
				alignment.print_result(100)
			i += 1
			print('-----------------------\n')
		print('--------------------Alignment END--------------------\n')
		again = input('Use the same seqs to align again Enter 0, Different sequence Enter 1, Exit Enter any key: ')
		if again == '0':
			go_align(lifex , ref)
		elif again == '1':
			lifex , ref = get_seqs()
			go_align(lifex , ref)
		else:
			a = input('Press any key to end the process')	
					
	except Exception as e:
		print(e)

if __name__ == '__main__':
	lifex , ref = get_seqs()
	go_align(lifex , ref )

