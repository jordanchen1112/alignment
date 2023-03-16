import pandas as pd
import random

def random_genertor(n):
	list1 = ['A','B','C','D','E','F','G','H']
	string = ''
	for i in range(n):
		b=random.randint(0,len(list1)-1)
		string += list1[b]

	return string

class Alignment:
	""" Alignment parameters """
	# i (seq) -> life X
	# j (seq_query)
	
	def __init__(self, seq: str, seq_query: str, gap_penalty , extend_gap_penalty, match_score, mismatch_score):
		self.seq = seq
		self.seq_query = seq_query
		self.gap_penalty = gap_penalty
		self.extend_gap_penalty = extend_gap_penalty
		self.match_score = match_score
		self.mismatch_score = mismatch_score
		self.align_start = [len(self.seq_query),len(self.seq)]
	
	def AA_matrix(self):
		# 胺基酸對應分數矩陣
		amino_acid = ['A','B','C','D','E','F','G','H']
		aa_matrix ={}
		for a in amino_acid:
			d = {}
			for a2 in amino_acid:
				if a2 == a:
					d[a2] = self.match_score  # 相同胺基酸分數
				else:  
					d[a2] = self.mismatch_score # 相異胺基酸扣分
			aa_matrix[a] = d

		return aa_matrix

	def create_score_matrix(self, affine): 
		# Create scoring matrix filled with 0 (int) to record the score
		score_matrix = []
		if affine == True:
			for j in range(len(self.seq_query) + 1):
				tmp =[]
				for i in range(len(self.seq) + 1):
					if i == 0 and j == 0:						
						tmp.append(0)
					elif i == 0 or j == 0:
						tmp.append(-100000)
					else:
						tmp.append(0)				
				score_matrix.append(tmp)
		elif affine == False:
			for j in range(len(self.seq_query) + 1):
				tmp =[]
				for i in range(len(self.seq) + 1):
					if i == 0 and j == 0:						
						tmp.append(0)
					elif i == 0:
						tmp.append(j * self.gap_penalty)
					elif j == 0:
						tmp.append(i * self.gap_penalty)
					else:
						tmp.append(0)				
				score_matrix.append(tmp)

		return score_matrix

	def create_xy_matrix(self):
		y_matrix = []
		x_matrix = []
		for j in range(len(self.seq_query) + 1):
			tmpx = []
			tmpy = []
			for i in range(len(self.seq) + 1):
				if i == 0 and j == 0:
					tmpx.append(0)
					tmpy.append(0)
				elif j == 0:						
					tmpx.append(-100000)
					tmpy.append(self.gap_penalty + self.extend_gap_penalty * (i-1))
				elif i == 0: 
					tmpy.append(-100000)
					tmpx.append(self.gap_penalty + self.extend_gap_penalty * (j-1))
				else:
					tmpx.append(0)
					tmpy.append(0)
			x_matrix.append(tmpx)
			y_matrix.append(tmpy)

		return x_matrix, y_matrix

	def create_trace_matrix(self):	
		# Create trace matrix filled with [0] (list) to record the arrow ; 
		# Score the first row and the first column

		trace_matrix = []
		for j in range(len(self.seq_query) + 1):
			tmp =[]
			for i in range(len(self.seq) + 1):
				tmp.append([0])
			trace_matrix.append(tmp)

		for i in range(len(self.seq) + 1):
			trace_matrix[0][i] = [0]
		for j in range(len(self.seq_query) + 1):
			trace_matrix[j][0] = [1]

		return trace_matrix	

	def score(self):
		pass

	def align(self):
		pass

	def find_next_step(self, matrix, trace_matrix,j,i):
		# 根據trace matrix 以及scoring matrix 回推位置[j,i]的胺基酸align情形
		# next_step = 0 表示從左至右(seq/-); next_step = 1表示從上至下(-/seq+query); 
		# next_step = 2 表往右下(seq/seq_query)

		possible_ways = {}
		next_step = 0
		back_ways = [matrix[j-1][i], matrix[j][i-1], matrix[j-1][i-1]]
		for ways in trace_matrix[j][i]:
			possible_ways[ways] = back_ways[ways]
			next_step = max(possible_ways.items(),key= lambda x:x[1])[0]
		return next_step

	def print_result(self, n = 40):
		result = self.result
		print(f'Align from {self.head} to {self.align_start}')
		headj = self.head[0] # seq_query
		headi = self.head[1] # seq

		if n > len(result[0]):
			n = len(result[0])
		# Interpret the alignment result
		for i in range(0,len(result[0]),n):
			n1 = len(result[0][i:i+n]) - result[0][i:i+n].count('-')
			n2 = len(result[2][i:i+n]) - result[2][i:i+n].count('-')

			print(f'\n{headj}   {result[0][i:i+n]}   {headj+n1}')
			space0 = ' ' * len(str(headj))
			space1 = ' ' * len(str(headj+n))
			print(f'{space0}   {result[1][i:i+n]}   {space1}')
			print(f'{headi}   {result[2][i:i+n]}   {headi+n2}\n')
			i += n
			headj += n1
			headi += n2
		return 0

class Global_Alignment(Alignment):

	def __init__(self, seq: str, seq_query: str, gap_penalty, extend_gap_penalty, match_score, mismatch_score):
		super().__init__(seq, seq_query, gap_penalty, extend_gap_penalty, match_score, mismatch_score)
		self.type = 'Global'
		self.align_start = [len(self.seq_query),len(self.seq)]

	def Score(self, affine, bounding):
		#Score the matrix from (0,0) to (m,n); record the score/arrow on the score/trace matrix
		# 建立score 和 trace matrix
		matrix = self.create_score_matrix(affine = affine)
		trace_matrix = self.create_trace_matrix()
		x_matrix, y_matrix = self.create_xy_matrix()
		aa_matrix = self.AA_matrix()

		all_max = -10000
		maxij = [len(self.seq_query),len(self.seq)]
		trace_back_start = [0,0]
		#填入score matrix的第一行和第一列
		for j in range(0, len(matrix), 1):			
			for i in range(0, len(matrix[0]), 1):
				if bounding:
					if i >= bounding * len(self.seq) + j:
						continue
					if j >= bounding * len(self.seq_query) + i:
						continue
	
				if i == 0 or j == 0:
					continue
				"""從(1,1)填入到(i,j)"""

				# 判斷胺基酸分數(aa_score)											
				aa_score = aa_matrix[self.seq_query[j-1]][self.seq[i-1]]

				# 比較三個方向的score大小(max_score:int), 並記錄trace matrix
				# trace back方向可能不只一個(0:'左' ; 1: '上' ; 2: '左上')

				# Affine Gap_Penalty
				if affine == True:			
					x_matrix[j][i] = round(max(x_matrix[j][i-1] + self.extend_gap_penalty, matrix[j][i-1] + self.gap_penalty, y_matrix[j][i-1] + self.gap_penalty), 2)
					y_matrix[j][i] = round(max(y_matrix[j-1][i] + self.extend_gap_penalty, matrix[j-1][i] + self.gap_penalty, x_matrix[j-1][i] + self.gap_penalty), 2)
					score_list = [x_matrix[j-1][i-1] + aa_score, y_matrix[j-1][i-1] + aa_score, matrix[j-1][i-1] + aa_score]
					matrix[j][i] = round(max(score_list), 2)
					max_score = max(x_matrix[j][i], y_matrix[j][i], matrix[j][i])
				else:
					score_list = [matrix[j][i-1] + self.gap_penalty, matrix[j-1][i] + self.gap_penalty, matrix[j-1][i-1] + aa_score]
					max_score = max(score_list)
					matrix[j][i] = round(max_score , 2)
					trace_matrix[j][i] = [i for i, x in enumerate(score_list) if x == max_score]

				#紀錄最大值及位置
				if i == len(self.seq) or j == len(self.seq_query):
					all_max = max_score if max_score >= all_max else all_max
					maxij = [j,i] if max_score >= all_max else maxij

		return matrix, trace_matrix, x_matrix, y_matrix, maxij

	def align(self, affine = False, bounding = None, stop_conditon = '&', align_start = 'End'):	
		# Align two sequences and return the result
		matrix, trace_matrix, x_matrix, y_matrix, maxij = self.Score(affine = affine, bounding = bounding)

		# 評分標準
		opengap = False
		score = 0
		the_same_score = 10/5
		dif_score = -8/5
		firstgap = self.gap_penalty / 5 
		extendgap = self.extend_gap_penalty /5 if self.extend_gap_penalty else self.gap_penalty

		# 紀錄align結果
		seq_query_record = ''
		seq_record = ''
		middle = ''

		# 決定Align開始位置(j,i)
		if align_start == 'End':
			start = [len(self.seq_query), len(self.seq)]
		elif align_start == 'max':
			start = maxij
		else:
			print('Error: Align start')
		self.align_start = start
		j = start[0]
		i = start[1]					
		print(f'Alignment type = {self.type}')	

		# 開始Align
		# 先決定下一步走哪個方向(find_next_step()),再填入align序列,最後才走到下一格
		# 利用->middle:str建立中間的align表示符號(相同:'|', 空格:'-')
		# 利用opengap = True/False 判斷屬於new gap或extend gap, 給予相異之分數
		while i >= 0 and j >= 0:
			if stop_conditon == '&':
				stop = i == 0 and j == 0
			elif stop_conditon == 'or':
				stop = i == 0 or j == 0
				# print(j,i,stop)
			elif stop_conditon == '0':
				s = max([x_matrix[j][i], y_matrix[j][i], matrix[j][i]]) if affine == True else matrix[j][i]
				stop =  s == 0
			if stop:
				break

			if affine == True and i == 0:
				next_step = 1
			elif affine == True and j == 0:
				next_step = 0
			elif affine == True:
				s = [x_matrix[j][i], y_matrix[j][i], matrix[j][i]]
				next_step = s.index(max(s))
			else:
				next_step = self.find_next_step(matrix, trace_matrix,j,i)

			if next_step == 2:
				if opengap == True:
					opengap = False
				seq_query_record = self.seq_query[j-1] + seq_query_record
				seq_record = self.seq[i-1] + seq_record
				if self.seq_query[j-1] == self.seq[i-1]:
					middle = "|" + middle
					score += the_same_score
				else:
					middle = " " + middle
					score += dif_score
				i = i - 1
				j = j - 1
			else:
				if opengap == False:
					opengap = True
					score += firstgap
				elif opengap == True:
					score += extendgap

				if next_step == 0: 			
					seq_query_record = '-' + seq_query_record
					seq_record = self.seq[i-1] + seq_record
					middle= " " + middle
					i = i - 1

				elif next_step == 1:			
					seq_query_record  = self.seq_query[j-1] + seq_query_record
					seq_record = '-' + seq_record 
					middle = " " + middle
					j = j - 1
			# matrix[j][i] = 'Go' + str(matrix[j][i])
		self.head = [j,i]
		result = [seq_query_record ,middle, seq_record]
		identity= round(list(result[1]).count('|') / len(result[0]) * 100, 2)
		self.result = result
		self.identity = identity
		self.score = score

		# Print similarity identity(%) and Alignment score  
		print(f'Alignnent similarity: {self.identity}%')
		print(f'Alignnent score: {self.score}')
		# print(pd.DataFrame(matrix))
		# print(pd.DataFrame(x_matrix))
		# print(pd.DataFrame(y_matrix))
		# print(pd.DataFrame(trace_matrix))
		return self.result

class SemiGlobal_Alignment(Global_Alignment):

	def __init__(self, seq: str, seq_query: str, gap_penalty, extend_gap_penalty, match_score, mismatch_score):
		super().__init__(seq, seq_query, gap_penalty, extend_gap_penalty, match_score, mismatch_score)
		self.type = 'Semi-Global'

	def Score(self, affine, bounding):
		return super().Score(affine = affine, bounding = bounding)

	def align(self, affine = False, bounding = None, stop_conditon = 'or', align_start = 'max'):
		super().align(affine = affine, bounding = bounding, stop_conditon = stop_conditon, align_start = align_start)
		return self.result

class Local_Alignment(Global_Alignment):

	def __init__(self, seq: str, seq_query: str, gap_penalty, extend_gap_penalty, match_score, mismatch_score):
		super().__init__(seq, seq_query, gap_penalty, extend_gap_penalty, match_score, mismatch_score)
		self.type = 'Local'

	def Score(self, affine, bounding):
		return super().Score(affine = affine, bounding = bounding)

	def align(self, affine = False, bounding = True, stop_conditon = '0', align_start = 'max'):
		super().align(affine = affine, bounding = bounding, stop_conditon = stop_conditon, align_start = align_start)
		return self.result

def main():
	seq = 'CBGECHEFDEBHGCFHDBFCBEAEBGHDBFCDADAHFEEACFAEBBEAEHACFDEABCFACFDHFBGHADHFFHGAEBDFGAAHHFGFCEABAABGHEBEABFFDEBHGAHGAHAFGABDABAGFBEAFCABAGFBEAFCEACBDFEAHGAGHDDBEBCFEAFBEFCBBEEABDAFHHEEABGEABGDHFCDABCFAGBCAGCHFAFGFBFBGHDBEGABAGHEBCEFEABGADHGECHGFCDBBEEAGBFHABDAEABEADHEEBHFGABDAEAFACBHEHHBDBEAFBEFCBBEAACFDFCGFEABGAEAFABGEFHFFFGEAHGAEAFAGBFHFDADABGFDFAECHFBFHDBFCBEAEBGHDBFCDAAGFEFDDBBABBBBBEAGHEAGFBEAFCDABDAAHDEFEBEFEABGABAGFFCFCAHGACFDFGEAABHAFACHGBEFABCEBDEFDAAACFFABBBBDACBAFEABEGACAABDAGHCFEEAFEABEGACAACBGAEABDABDABAAEBFDBCEFABGEFCACFEBEBHGAHBGFGAEAFAHFBEEAAHGAEBEBAEABEAABDABDDFFFEBEFEAHGFCAEAFAABDEACFABFBCDAHGAEAFADHFAEFAABCDABEFDEFCFAHGADHEEBHFGHFDAGBCFCDABGAEAFAEFCFBEABGEADFCEFCFBEAEBBFCDAHGAEAFADDBGABGABAHBEFACBGHFAHGAGFCEFCCBEFDEAFGHGEADABCDDAAFHEEBABBHHDAHBBGHCBHAEAFEABEGABBHABEACHGBAGBDAAAAFCCBGDABGEAAFCCBGDABBAGBEAEHEAABGDAAABCDEABBBGBEADGBDFDAACBBGFABBAABEAEHHHFCAFBEAEFCEEFDAAAFCDHGBEAHCDFCGBEBHGBEAFCHEFEFDAAGCHEBDAABGEADDAFBEABBBBBEABGEABGFCBGDAAHCFGFGAFEABEGABBBFBGABGEFHFFFGEBEADECFDEFCFDABGAFBFFHEADAADFDABCFGBABGEABEFCBHBFHGABBABBABGEABDAEABHDBFCDAAEBGHABFAADBEDDEAACHGBEBGGABBHEAFCAHBBAEHAGBHAEABHBBGDEAECFHCCFDDDEBGDFGAACCAAGFCEAFCEAEAFBADFBEEBAADBAAAFAEABDDGAAHHFAABBEFDEAACHGBEBGHABBHEAFCAHBBAEHGGBHEEABHBBGDEAECFHACFDBDEBGDFGAAACAAGFCEAFEEAEAFBABEFEBAACFDFEFAAHEFGEBBEAGHDDBEBCBEBHGAGHEEHHBGHAEAFAEFBEAAHGABGABGBFBEGC'
	seq_query = 'CAGAAEADEFEFGECAHBGHEADAFFDABGHADADBCDEFEFGEABECABHAAACBACEBEFCAGHGFFCFCACFEACABBCDAFFDAAEBEEFCAEBBHBGAHABEFEAEHEAABGAFDHEHHBDBEAACFDFCGBEBHGACDAFBDFCCAECGADAHFEAEBFGFDBBGHCCECGADAHFEAEBFGFDBBGHAABDADFCGBDFEABGABGDEBEFEFAHGAFDHEHHBABGEAFGHEFEBHGBCBAHGAEBBHBGAGBEBHGBEAFGBGFCDBEBAGHCAFBGBABFBCDGAAHGAGHGFFCFCACFEACABBEADAFAFBEFABAGFCBABGEFCFDEBGHADAFFDAGAECGADAHFAFBCHCAEFDDCBCFEAEAFAEBDECBCFEBHGEAABCBEBEEACFABGBHCEAFDHEHHBDBEAACFDFCGBEBHGAHGADABGFDFAHABEFAEHEAABGGAGCHFAAFCADAFFDAEABAHHEAFBGBAAHBGEDABCHFEAEBBHBGAHABEFEAEHEAABGAFDHEHHBABGEAACFDFCGBEBHGGAAFCFABDAFBACFGEFDEBHGAHGAEABDADAFFDAGCCDDBFGEBGBDAGBFFAHGADABGFDFAHABEFAEHEAABGABDADHFDBADABGFGDBDEAHABDAAEBDECBCFEFABGADHBDEBEAHBEFCDAHGADABGFDFEAEBBHBGEAEAFADHFEAFBDEAABDBGBDAHDFBGADFCCHFGEBGHADHFGECBFDEABFDECBBEADHFEAABGCBDBEABGEABGEBBGGAAEAFBAEHGBEAEBDFACBHAHABEFDADBGAECBGFEAEHGHAHBBAEHAEBGGFCFGEAABCBEDGAEAFBACFDEAEBGFAEAFADHBDEBEADABEEHHAHBEFCDGAAHAFGAEAFBABCFAEFGABFBCDAHEEAHABDAAABGFABCBEBEBAHGAACHEFDEBHGAGAAEBDFAHEAFCAHABEFDABGEAEHEAABGDAFDFAEAFACBHFDHGBCAEHAEFGBGFAEAFAAHDBEBHGABGEADHFFFGBDBEFAHBEAAHEAFCDEADABGFDFAHABEFAEHEAABGADBGADHFFFGBDBEFAHBEAAHEAFCDACBAHABDEEFACFHBDEFCEABADAFDBBEADHFGEABGAEAFAHBEFCGACAAACEBBHBGAABDADAFDBBEAHFHHCBAABDBEAFGGBCHGFFGEEAEAFAHFDEADHBDEABDADABEEHHAHBEFCABGEAFBDEABDAEFFAAHBEFCEAEAFDADABGFDFAHABEFAEHEAABGACFDEAEBGFABGAHFDEADHBDEBEAHBEFCDGAAGCHFACAAHAEHACABAEAECGADAHFAHCHFAAGHFGEABCHFEAHAAEHAAFADABGFDFAHABEFAEHEAABGDABGAGBFEEABGGFDEBHBEBHGAHGABGAEAFAHFDEADHBDEAHGAEBBHBGGACCAAAAFFBGAHGFCGBDABGHAEFDCFBDFAEAFAGHHEEAAFECHDAFFBDBEABGEFDECBAEBDDABCHFAEAFAAFBGBAFFGEBEAAAHEAHCHBGBDFAHEADFEBAEHAFBDFACBHBDDFFFEBEBHGDAADAFGEAFFGHFADBFGBACAACBEAEAFAHBEFCAAHEEFEBHGAEFDECHBAEAFAABCBEDEAEFEABEAGBCFDABGEAABEAHHFGBFDEABEEAHGAEAFDFAGBDEHCDAFBDFADABGFDFAHABEFAEHEAABGDAEHACFDHFFABGAFGEBGHFCFEADAFDBFDGACCAAAADABGFDFAHABEFAEHEAABGAFDHEHHBDBEAACFDFCGBEBHGABDADCBEBDBEABGAEAFAHHCEEGAEBBHBGAHHGFCGFFGEADAHFEEAHGGFCAEAFAACFDFCGBEBHGAAHEBDBEABGDEFEFAACHABCBEAEAFAFFBDDBHGAHGAABCFGFEAHBDEFAEHAEAFAHBEFCEABGHBEAEHACFBEEAEAFAAFECHDAFFBDBEABGEFDECBAGFBCAHFDEADHBDEABGEAHGFCGBDABGHEAEABEAHBEEAACHEFDEAEAFDFAFGEBGHFCFEADAFDBFDCBGEFCADAFFDAEABACFBFFDEABABFFDEBHGAEABEABDAEABDAHFADGHHAAFFBGAABGFAGBCBHFDAEHDBEAEBGHFBHFDABGAEBGGFCFGEAAEBDFDAHGAEAFAFBCEAGAEHAEAFADABGFDFAHABEFAEHEAABGDAABGFAEAFADBFFACFABGBHCABDAAFFBGHFCCAAAAAECGADAHFADBBEEABDEFBEEBEAEABDABDAGFCBABGEFCFDEBGHABFFDEBHGGADHAGBCABDABADGHHEACDABGFDFAHABEFAEHEAABGAFBBCFAEBDFAAFFBGAHAHAABGFAEBGGFCFGEAEBGHFBHFDGACGCHFAEAFADFCCFGEADFBCDAAHGABFDECBBEADABGFDFAHABEFAEHEAABGDAABGFAEAFAEHGHFCAHBGFEFGHEAAHGAHABDEEFACFHBDEFCAEABGAHEAFCAAEBDFDGAEABEAFFBGDAHFHHCBAABDBEABDHEBEBHGAFBDFDAEBGGFCFGEAFGHEFEBHGGAEBGHFBHFAEFGFEHAFFGEABDAEBDFAEABEGCCGFAEADFFFDEFCADHFCDFADFHHFDEBHGCBGAEABDADFFFDEFCEABAHHEAFBGBAAHBGEDAGCHFAEAFDFAGBDFAEBEDDGAHBGBGHAEAFADHGDFCGFEABGGHCFBEBHGAHCAABAFCDACFGHCFAEAFADAFFDAFDAEABEAHBEEABGDCFBDFAEAFAFGEFCDEBGEBGHABGEABCDHCAEBHGAHGADEFEFGEDGACFDBEFDEABGAEAFDFADAFFDAFDADBGACFACFDHCEABGEAAFEABGEHAEAFAEFBCGBGHAAEBEGHCFAEABEAHBEEAACHFHEFADEFEFGEDABGEFCADEBDDFDGCAAACCCFGFCFGDFCDAFGEAFGFAGEADGFDGADABAEADGAEGADAHFABGEAEGFDGADAHFACAACGFFCDFCBECHCHBGBDFFFCDFCBABGEADFEFGBFFABGADFBEEADFEBDFBGDABGAEBBHBGFDFAHBEFCDGAFBCBGFAAHEEFEBHGACFEEFEBGAEFCCDHFCEFGC'

	gap_penalty = -4
	extend_gap_penalty = gap_penalty/10
	match_score = 10
	mismatch_score = -8

	alignment = Local_Alignment(seq, seq_query, gap_penalty,extend_gap_penalty, match_score, mismatch_score)
	alignment.align(affine = True)
	alignment.print_result()


if __name__ == '__main__':
	main()


