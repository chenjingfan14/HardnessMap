import os, shutil

myFiles = []

NewPath=os.path.normpath("C:/FullFileTargetPath")

BasePath=os.path.normpath("C:/FullExistingTargetPath")

# iterate through all the lines
for line in open('MyTargetResultsFile.spe').readlines():
    # skip this line if it doesn't look like an assignment
	if not BasePath in line: continue
	
	left, right = line.split('>', 1)
	# print right
	left, right = right.split('<',1)
	print left
	if '.jpg' in left:
		myFiles.append(left)

for entry in myFiles:
	#Copy the other photos, change trailing file names b/c of split op
	#File name ends in '_Result.jpg', change to '_Result_Thumb.jpg' to each
	#File name ends in 'x_Result.jpg', need to change to 'x.jpg'
	thumb = entry.replace('_Result.jpg','_Result_Thumb.jpg')
	raw = entry.replace('x_Result.jpg','x.jpg')
	#copy the file to another directory
	try:
		shutil.copy2(entry, NewPath)
		shutil.copy2(thumb, NewPath)
		shutil.copy2(raw, NewPath)
	except IOError, e:
		if e.errno == 2:
			#doesn't exist
			pass
		else:
			#uhoh
			raise

