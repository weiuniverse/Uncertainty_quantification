import tsv
'''
truth = tsv.TsvReader(open("../data/poly_mo/poly_truth.tsv"))
count1 = 0
for parts in truth:
  # Here parts is a list of strings, one per tab-separated column.
  # Make sure you handle not having enough fields, or not being able to
  # parse numbers where you expect them.
  # print("Record with {} fields: {}".format(len(parts), parts))
  count1 = count1 + 1
  # print parts

quant_bootstraps = tsv.TsvReader(open("../data/poly_mo/quant_bootstraps.tsv"))
count = 0
for parts in quant_bootstraps:
  # Here parts is a list of strings, one per tab-separated column.
  # Make sure you handle not having enough fields, or not being able to
  # parse numbers where you expect them.
  # print("Record with {} fields: {}".format(len(parts), parts))
  count = count + 1
  # print len(parts)
print count
'''
quant = open("../data/poly_mo/quant.sf")
lines = quant.readlines()
# print l
for line in lines:
    line = line[:-1]
    l = line.split('\t')
    print len(l)
    if len(l) <> 5:
        print "error"
        break
# quant  = tsv.TsvReader(open("../data/poly_mo/quant.sf"))
# for parts in quant:
#     print quant