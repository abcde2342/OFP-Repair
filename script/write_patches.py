import io, os, struct

def doubleToInt(f):
    return struct.unpack('<Q', struct.pack('<d', f))[0]
def intToDouble(i):
    return struct.unpack('<d', struct.pack('<Q', i))[0]

def fmtFP(x):
	return "{:.17e}".format(x)

class PatchInfo:
	data_path   = './data/polynomial.dat'
	funcID_path = './data/funcID.dat'
	valid_path  = './data/validation.dat'

	def __init__(self):
		self.check_file(PatchInfo.data_path)
		self.check_file(PatchInfo.funcID_path)
		self.read_data()

	def check_file(self, filepath):
		if not os.path.exists(filepath):
			print('Missing data file:', filepath, ". Run [bin/gsl-fit.out] to generate it.")
			exit()
		return True

	def read_data(self):
		with io.open(PatchInfo.data_path) as f:
			self.rawdata = f.readlines()
		with io.open(PatchInfo.funcID_path) as f:
			self.patch_no = f.readline().split(' ')[-1].strip()
			# self.patch_no = int(f.readline()) + 1 # start from 1
		with io.open(PatchInfo.valid_path) as f:
			self.patch_valid = int(f.readline())
		self.polymode = int(self.rawdata[0], base=16)
		assert (self.polymode==1 or self.polymode==2), "Template for power series not found."
		i64lb = int(self.rawdata[1], base=16)
		i64rb = int(self.rawdata[2], base=16)
		self.patch_lb = intToDouble(i64lb)
		self.patch_rb = intToDouble(i64rb)
		i64c1 = int(self.rawdata[3], base=16)
		i64c2 = int(self.rawdata[4], base=16)
		self.patch_c1 = intToDouble(i64c1)
		self.patch_c2 = intToDouble(i64c2)

		self.patch_coefs = []
		for line in self.rawdata[5:]:
			i64x = int(line, base=16)
			fpx  = intToDouble(i64x)
			self.patch_coefs.append(fpx)
		self.patch_terms = len(self.patch_coefs)

	def print_patch_info(self):
		print('No.', self.patch_no)
		print('Valid?', self.patch_valid)
		print('lbound', fmtFP(self.patch_lb))
		print('rbound', fmtFP(self.patch_rb))
		print('terms ', self.patch_terms)
		print('coefs ')
		for coef in self.patch_coefs:
			print(fmtFP(coef))

class TemplateInfo:
	template_path = './inputs/cheby.template'
	patch_path    = './patches/patches.h'

	def __init__(self):
		self.read_template()
		if not os.path.exists('./patches'):
			os.mkdir('./patches')

	def read_template(self):
		with io.open(TemplateInfo.template_path) as f:
			lines = f.readlines()
			for index, line in enumerate(lines):
				if '__TEMPLATE_START__' in line:
					self.template_start_no = index
					self.binding_part = lines[0:index]
				if '__TEMPLATE_END__' in line:
					self.template_end_no = index
					self.template_part = lines[self.template_start_no+1:self.template_end_no]

	def write_binding_part(self):
		if not os.path.exists(TemplateInfo.patch_path):
			with io.open(TemplateInfo.patch_path, 'w') as f:
				for line in self.binding_part:
					f.write(line)

	def write_template_part(self, patchinfo):
		newlines = []
		for line in self.template_part:
			if '__PATCH_NO__' in line:
				line = line.replace('__PATCH_NO__', str(patchinfo.patch_no))
			if '__PATCH_IS_VALID__' in line:
				if patchinfo.patch_valid == 1:
					comment_str = 'Valid Patch'
				else:
					comment_str = 'Invalid Patch'
				line = line.replace('__PATCH_IS_VALID__', comment_str)
			if '__PATCH_LEFT_BOUND__' in line:
				line = line.replace('__PATCH_LEFT_BOUND__', fmtFP(patchinfo.patch_lb))
			if '__PATCH_RIGHT_BOUND__' in line:
				line = line.replace('__PATCH_RIGHT_BOUND__', fmtFP(patchinfo.patch_rb))
			if '__PATCH_TERMS__' in line:
				line = line.replace('__PATCH_TERMS__', str(patchinfo.patch_terms))
			if '__PATCH_COEFFICIENT__' in line:
				line = ''
				for coef in patchinfo.patch_coefs:
					line += '    '+fmtFP(coef)+',\n'
			if '__PATCH_ACC_MODE__' in line:
				line = line.strip() + " " + str(patchinfo.polymode) + '\n'
			if '__PATCH_ACC_BOOL__' in line:
				if patchinfo.polymode == 2:
					line = line.replace('__PATCH_ACC_BOOL__', 'true')
				else:
					line = line.replace('__PATCH_ACC_BOOL__', 'false')
			if '__PATCH_ACC_C1__' in line:
				if patchinfo.polymode == 2:
					line = line.replace('__PATCH_ACC_C1__', fmtFP(patchinfo.patch_c1))
				else:
					line = line.replace('__PATCH_ACC_C1__', '0.0')
			if '__PATCH_ACC_C2__' in line:
				if patchinfo.polymode == 2:
					line = line.replace('__PATCH_ACC_C2__', fmtFP(patchinfo.patch_c2))
				else:
					line = line.replace('__PATCH_ACC_C2__', '0.0')
			
			newlines.append(line)
		
		with io.open(TemplateInfo.patch_path, 'a') as f:
			for line in newlines:
				f.write(line)

	def write_patch(self, patchinfo):
		self.write_binding_part()
		self.write_template_part(patchinfo)

def main():
	pi = PatchInfo()
	
	ti = TemplateInfo()
	ti.write_patch(pi)

if __name__ == '__main__':
	main()