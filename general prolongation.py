from sympy import *
import prolong_Fast

# This is tedious and can be automated using itertools, but we are essentially initializing our variables
x,t,u,ux,ut,uxx,uxt,utt,utx,uxxx,uxtx,utxx,uxxt,uxtt,utxt,uttx,uttt,xi,xit,xiu,xix,xixx,xixt,xitx,xiux,xixu,xiuu,phi,phit,phiu,phix,phixx,phixt,phitx,phiux,phixu,phiuu,tau,taut,tauu,taux,tauxx,tautx,tauxt,tauux,tauxu,tauuu\
= symbols("x t u ux ut uxx uxt utt utx uxxx uxtx utxx uxxt uxtt utxt uttx uttt xi xit xiu xix xixx xixt xitx xiux xixu xiuu phi phit phiu phix phixx phixt phitx phiux phixu phiuu tau taut tauu taux tauxx tautx tauxt tauux tauxu tauuu")

# Here we define two polynomial rings that we will work over
R = PolynomialRing(ZZ,'x,t,u,ux,ut,uxx,uxt,utt,utx,uxxx,uxtx,utxx,uxxt,uxtt,utxt,uttx,uttt')
M = PolynomialRing(ZZ,'phit,phiu,phiuu,phixu,phixx,taut,tauu,tauuu,taux,tauxu,tauxx,xit,xiu,xiuu,xix,xixu,xixx')

f = Poly(phi-xi*ux-tau*ut, domain = R)

# Here we are defining our dependent variables and their jet coordinates
U = [u,ux,ut,uxx,uxt,utt,utx,uxxx,uxtx,utxx,uxxt,uxtt,utxt,uttx,uttt]

# Here we handle the order of mixed derivatives being the same
Subs = {utx:uxt, uxtx:uxxt, utxx:uxxt, utxt:uxtt, uttx:uxtt,xitx:xixt, xiux:xixu, phitx:phixt, phiux:phixu, tautx:tauxt,tauux:tauxu}

# This is a simple symbolic differentiation method
def d(expr, Var):
    if expr == 1:
        return 0
    else:
        return Poly(str(expr)+str(Var)+"+"+str(expr)+"u"+"*"+"u"+str(Var), domain = R)

# U contains all dependent variables and their derivatives that show up in the polynomial.
# This is used in substituting expresions like Derivative(u(x),x) and u(x) with ux or u respecively.
def Der(poly, Var, U):
    g = 0
    n = len(poly.coeffs())
    for i in range(n):
        g = g + d(LM(poly),Var)*LC(poly)+LM(poly)*LC(poly).subs({u:u(Var) for u in U}).diff(Var)
        poly = poly-LC(poly)*LM(poly)
    g = g.subs({Derivative(u(Var), Var): str(u)+str(Var) for u in U})
    g = g.subs({str(u)+"("+str(Var)+")": str(u) for u in U})
    g = g.subs(Subs)
    return g

if __name__ == "__main__":
    
    # Example 1
    # g is phi^xx and h is phi^t
    g = expand(Der(Poly(Der(Poly(phi-xi*ux-tau*ut,domain=R),x,U),domain=R),x,U)+xi*uxxx+tau*uxxt)
    h = expand(Der(Poly(phi-xi*ux-tau*ut,domain=R),t,U)+xi*uxt+tau*utt)
    # f is the difference of g and h, and we substitute ut with uxx (System is satisfied)
    f = g-h
    f = f.subs({ut:uxx})
    # Here we convert f into a polynomial in the ring M
    r = Poly(f,domain = M)
    print "Here is the system of odes for determining the symmetries of the heat equation:\n"
    print '\n'.join(str(p)+" = 0" for p in r.coeffs()) 