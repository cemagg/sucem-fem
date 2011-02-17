import csv

def read_multi(iter, no):
    while True:
        multi = []
        for i in range(no):
            try:
                multi.append(iter.next())
            except StopIteration:
                if i > 0 : yield multi
                raise StopIteration
        yield multi
        

def read_grad_lambdas(fileh):
    csvreader = csv.reader(fileh)
    grad_lambdas = []
    for el in read_multi(csvreader,4):
        grad_lambdas.append([
            [float(comp) for comp in lamb]
            for lamb in el
            ])
    return grad_lambdas
                        
