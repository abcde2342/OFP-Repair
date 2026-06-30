import io

tmpfiles = [
    '/tmp/1.log',
    '/tmp/2.log',
    '/tmp/3.log',
    '/tmp/4.log',
    '/tmp/5.log',
    '/tmp/6.log',
    '/tmp/7.log',
    '/tmp/8.log',
    '/tmp/9.log',
    '/tmp/10.log'
]

count = 0
for f in tmpfiles:
    print('Runtime cost of Aceso on round', count+1)
    count += 1
    with io.open(f) as ff:
        data = ff.readlines()
    lines = []
    for line in data:
        if line.startswith("Time on Synthesizing"):
            lines.append(line)
    for line in lines:
        print(line.split()[-2])
    print("=================\n\n")