


import sys

def main():
    comp_ops = {'>' : (lambda x,y: x >  y), 
               '>=' : (lambda x,y: x >= y),
               '<'  : (lambda x,y: x <  y),
               '<=' : (lambda x,y: x <= y),
               '=='  : (lambda x,y: x == y),
               '!=' : (lambda x,y: x != y)}
    cmp_type = sys.argv[1]
    if(cmp_type not in comp_ops):
        print( "operator {} is not defined. Use one of the following".format(cmp_type),file=sys.stderr)
        for i in comp_ops.keys():
            print(i,file=sys.stderr)
        exit(-1)
    required = sys.argv[2]
    toolvers = sys.argv[3]
    
   
    MAXFIELD = 10000
    
    refields = reversed(required.split("."))
    tofields = reversed(toolvers.split("."))

    reqvalue = 0

    for i,v in enumerate(refields):
        reqvalue += int(v) * (MAXFIELD ** i)

    toovalue = 0
    for i,v in enumerate(tofields):
        toovalue += int(v) * (MAXFIELD ** i)

    return 0 if ((comp_ops[cmp_type](toovalue,reqvalue))) else 1
if __name__ == '__main__':
    sys.exit(main())
