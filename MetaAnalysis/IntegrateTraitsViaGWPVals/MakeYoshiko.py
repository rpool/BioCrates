#!/usr/bin/python
# script mailed by Mohammed on 10/02/2014
import sys
f = open(sys.argv[1])
t = float(sys.argv[2])

# read labels
l = f.readline().rstrip("\n").split(",")[1:]

# read matrix
d = []
for line in f:
  d.append([])
  s = line.rstrip("\n").split(",")[1:]
  for val in s:
    d[-1].append(float(val))

print len(l)
for label in l:
  print label

for i in range(len(l)-1):
    print " ".join(map(lambda x: str(x - t), d[i][i+1:]))