o.alpha=m.alpha
o.beta=m.beta
o.prob=m.prob

terms_to_be_added=[85, 1935,679, 760, 1162]
seed=1

o.initialize()
m.initialize()

m.percentage_parameter_change=1
o.percentage_parameter_change=1

included=[]
for term in terms_to_be_added:
    print "m",m.n00,m.n01,m.n10,m.n11,m.get_score()
    print "o",o.n00,o.n01,o.n10,o.n11,o.get_score()

    if term in included:
        included.remove(term)
        o.remove_set(term)
        m.remove_set(term)
    elif term>=0:
        included.append(term)
        o.add_set(term)
        m.add_set(term)
    else:
        random.seed(seed)
        m.propose_state(100)
        random.seed(seed)
        o.propose_state(100)
        seed=random.random()*1000000

print "m",m.n00,m.n01,m.n10,m.n11,m.get_score()
print "o",o.n00,o.n01,o.n10,o.n11,o.get_score()

m.percentage_parameter_change=1
o.percentage_parameter_change=1

random.seed(seed)
m.propose_state(100)
random.seed(seed)
o.propose_state(100)

print "m",m.n00,m.n01,m.n10,m.n11,m.get_score()
print "o",o.n00,o.n01,o.n10,o.n11,o.get_score()
