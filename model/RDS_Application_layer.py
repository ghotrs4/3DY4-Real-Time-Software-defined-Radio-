def process_rds_data(msgs, prevPTYcode, prevPIcode, count):
	print ("A block", msgs.a)
	print ("B Block", msgs.b)
	print ("C Block", msgs.c)
	print ("D Block", msgs.d)
	PIcode = [] #program identification code in hex
	group_type = [] 
	PTYcode = '' #program type code
	DCBcode = [] #decoder controller bits
	PSNS = ['']*8 #program service name segment
	Program_Type_Codes = {
		'00000': 'No programme type or undefined',
		'00001': 'News',
		'00010': 'Current Affairs',
		'00011': 'Information',
		'00100': 'Sport',
		'00101': 'Education',
		'00110': 'Drama',
		'00111': 'Culture',
		'01000': 'Science',
		'01001': 'Varied',
		'01010': 'Pop Music',
		'01011': 'Rock Music',
		'01100': 'Easy Listening Music',
		'01101': 'Light classical',
		'01110': 'Serious classical',
		'01111': 'Other Music',
		'10000': 'Weather',
		'10001': 'Finance',
		'10010': "Children's programmes",
		'10011': 'Social Affairs',
		'10100': 'Religion',
		'10101': 'Phone In',
		'10110': 'Travel',
		'10111': 'Leisure',
		'11000': 'Jazz Music',
		'11001': 'Country Music',
		'11010': 'National Music',
		'11011': 'Oldies Music',
		'11100': 'Folk Music',
		'11101': 'Documentary',
		'11110': 'Alarm Test',
		'11111': 'Alarm'
	}
	programServiceName = {
		'0100 0001':'A',
		'0100 0010':'B',
		'0100 0011':'C',
		'0100 0100':'D',
		'0100 0101':'E',
		'0100 0110':'F',
		'0100 0111':'G',
		'0100 1000':'H',
		'0100 1001':'I',
		'0100 1010':'J',
		'0100 1011':'K',
		'0100 1100':'L',
		'0100 1101':'M',
		'0100 1110':'N',
		'0100 1111':'O',
		'0101 0000':'P',
		'0101 0001':'Q',
		'0101 0010':'R',
		'0101 0011':'S',
		'0101 0100':'T',
		'0101 0101':'U',
		'0101 0110':'V',
		'0101 0111':'W',
		'0101 1000':'X',
		'0101 1001':'Y',
		'0101 1010':'Z',

		'0110 0001':'a',
		'0110 0010':'b',
		'0110 0011':'c',
		'0110 0100':'d',
		'0110 0101':'e',
		'0110 0110':'f',
		'0110 0111':'g',
		'0110 1000':'h',
		'0110 1001':'i',
		'0110 1010':'j',
		'0110 1011':'k',
		'0110 1100':'l',
		'0110 1101':'m',
		'0110 1110':'n',
		'0110 1111':'o',
		'0111 0000':'p',
		'0111 0001':'q',
		'0111 0010':'r',
		'0111 0011':'s',
		'0111 0100':'t',
		'0111 0101':'u',
		'0111 0110':'v',
		'0111 0111':'w',
		'0111 1000':'x',
		'0111 1001':'y',
		'0111 1010':'z',

		'0010 1100':"'",
		'0010 1101':'-',
		'0010 1110':'.',
		'0010 1111':'/',

		'0011 0000':'0',
		'0011 0001':'1',
		'0011 0010':'2',
		'0011 0011':'3',
		'0011 0100':'4',
		'0011 0101':'5',
		'0011 0110':'6',
		'0011 0111':'7',
		'0011 1000':'8',
		'0011 1001':'9',
	}
	DCBLookup = {
		'00': 0,
		'01': 2,
		'10': 4,
		'11': 6
	}
	# Process the RDS data segments only if they are not empty
	if msgs.a:
		for i in range(0, 16, 4):
			binary_group = msgs.a[i:i+4]
			decimal_value = sum(binary_group[j] * (2 ** (3 - j)) for j in range(4))
			hex_digit = hex(decimal_value)[2:].upper()
			PIcode.append(''.join(hex_digit))
		PIcode = ''.join(PIcode)


	if msgs.b:
		group_type = msgs.b[0:5] #group type code
		#Dont need to extract TP bit
		PTYcode = ''.join(map(str, msgs.b[6:11])) #program type code
		#Dont need to extract TA and M/S bits
		DCBcode = msgs.b[14:16] #decoder controller bits
	
	if msgs.c:
		pass
	
	if msgs.d:
		group_parsed = ''.join(map(str, msgs.d[0:5]))
		if group_parsed == "00000":
			char1_coded = ''.join(map(str, msgs.d[0:8]))
			char2_coded = ''.join(map(str, msgs.d[8:16]))
			if char1_coded in programServiceName:
				char1 = programServiceName[char1_coded]
			else:
				char1 = ''
			if char2_coded in programServiceName:
				char2 = programServiceName[char2_coded]
			else:
				char2 = ''


			DCBcode = ''.join(map(str, msgs.d[0:2]))
			idx=DCBLookup[DCBcode]
			PSNS[idx] = char1
			PSNS[idx+1] = char2
			count +=1
			
		pass
			
	
	if (count == 4):
		Program_Service = ''.join(map(str, PSNS))
		print("Program service:", Program_Service) #print the program service
		count = 0

	#print(char1)
	if (PTYcode != prevPTYcode and PIcode != prevPIcode and PTYcode != '' and PIcode != []):
		if PTYcode in Program_Type_Codes:
			print("PI code:", PIcode)
			print("Program type:", Program_Type_Codes[PTYcode])
	
	return PTYcode, PIcode, count

if __name__ == "__main__":

	# do nothing when this module is launched on its own
	pass