data = open('doc-index-structure.txt')

pages = []
funcs = []

while 1:
	line = data.readline()
	if not line: break

	ds = line.split('\t')

	pages.append(ds[0])
	funcs.append(ds[1])
	
data.close()


out = file('doc-index.txt','w')

i = 0
for p in pages:

	f = funcs[i]

	if p == 'BTypes':
		line = '<li><codetype><a href="doc-files/'+p+'.html#'+f+'" title="'+p+'.jl">::'+f+'</a></codetype></li>\n'
	else:
		line = '<li><codefun><a href="doc-files/'+p+'.html#'+f+'" title="'+p+'.jl">'+f+'</a></codefun></li>\n'

	out.write(line)
	i += 1

out.close()