import shutil

g = open("submit", 'w')
h = open("combine", 'w')
h.write("cat ")
for i in range(10):
    val=str(2.5 + i*0.05)
    cutoff = 3.448*float(val)
    shutil.copyfile("in.water",(val+"/in.water"))
    shutil.copyfile("sub.sh",(val+"/sub.sh"))
    inp=val+"/pairstyle.inc"
    f = open(inp, 'w')

    f.close()

    g.write('cd %s\n' % val)
    g.write('sed -i -e "s@AAA@%s@g" in.water\n' % cutoff)
    g.write('msub sub.sh\n')
    g.write('cd ../\n')
    
    h.write('%s/henry.avg ' % val)
    
h.write(" > combine_henry")
h.close()
g.close()

