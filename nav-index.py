data = open('nav-structure.txt')

pages = []
funcs = []

while 1:
	line = data.readline()
	if not line: break

	ds = line.split('\t')

	pages.append(ds[0])
	funcs.append(ds[1])
	
data.close()


out = file('nav-index.txt','w')

i = 0
pp = "a"
for p in pages:
	if p != "*":
		if i != 0:
			out.write('</ul>\n</div>\n</li>')
			
		f = funcs[i]
		out.write('<li><a href=\\"'+p+'.html#\\">'+p+'.jl<span class=\\"caret\\"></span></a>\n<div>\n<ul>')
		out.write('<li><a href=\\"'+p+'.html#'+f+'\\">'+f+'</a></li>\n')
		pp = p

	else:
		f = funcs[i]
		line = '<li><a href=\\"'+pp+'.html#'+f+'\\">'+f+'</a></li>\n'
		out.write(line)

	i += 1

out.write('</ul>\n</div>\n')
out.close()
