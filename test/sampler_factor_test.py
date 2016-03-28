#This is the test code for sampler_factor


def test_cr(v, t, alpha, lambda_u):
	print "hello world"


def test_mw(v, t, alpha, lambda_u):
	#print "hello mengqing"
    d1 = v.shape
    d2 = t.shape
    for i in range(d1[0]):
        for j in range(d2[0]):
            vt_vector = np.multiply(v[i], t[j])
            precision_matrix = np.add(precision_matrix, alpha * np.dot(np.array([vt_vector]).T, np.array([vt_vector])))
    return precision_matrix


def test_main():
	print "The main function"