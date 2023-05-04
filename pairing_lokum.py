from sage.all_cmdline import *   # import sage library

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.misc.functional import cyclotomic_polynomial
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.schemes.elliptic_curves.ell_point import *

# this is much much faster with this statement:
# proof.arithmetic(False)
from sage.structure.proof.all import arithmetic

from pairing import *
from pairing_cp6_782 import *
from test_pairing import *
    
def check_curve_order(E, order):
    # randomized probabilistic check, usefull for large parameters
    ok = True
    i=0
    while ok and i < 10:
        P = E.random_element()
        ok = order*P == E(0)
        i += 1
    return ok

# Originally -> def test_lokum(u, a, b, q, r, tr, c, y, D, xi):
def test_lokum(a, b, q, r, tr, c, D):
    print("###################################################################")
    print("test curve E: y^2 = x^3 + a*x + b")
    print("a = {} b = {:#x}".format(a, b))

    #assert q == (tr**2 + D*y**2)//4
    print("q = {:#x} # {} bits".format(q, q.nbits()))
    Fq = GF(q)
    E = EllipticCurve([Fq(a),Fq(b)])
    P = E.random_point()
    print("P is a random_point on the elliptic curve!")
    print(P)
    print("The order of the curve is:", c)
    order_P_Q = r * c
    print("P * order_P_Q is")
    print(P * order_P_Q)
    res = check_curve_order(E,r*c)
    print("check_curve_order(E,r*c): {}".format(res))
    print("ate pairing, Miller loop scalar is tr-1 =\nT = {}\nT = {:#x}".format(tr-1, tr-1))

    ## THIS IS FOR EXTENDING LOKUM FIELDS 1-4-8
    xi = ZZ(2)
    Fqz = Fq['z']; (z,) = Fqz._first_ngens(1);
    Fq8_abs = Fq.extension(z**8 - xi, names=('j',)); (j,) = Fq8_abs._first_ngens(1)
    Fq4 = Fq.extension(z**4 - xi, names=('j',)); (j,) = Fq4._first_ngens(1);
    Fq4Z = Fq4['Z']; (Z,) = Fq4Z._first_ngens(1);
    Fq8 = Fq4.extension(Z**2 - j, names=('s',)); (s,) = Fq8._first_ngens(1)

    # EXAMPLE
    print("FQ8 is")
    F8 = GF(q**8, 'b')
    b = F8.gen()
    E8 = EllipticCurve(F8, [F8(a), F8(b)])
    print(E8)
    # EXAMPLE
    EK = E.base_extend(Fq8_abs)
    print(EK)
    
    print("A random point in the extension field, Q is: ")
    Q = EK.random_point()
    print(Q)

    print("A random point in the extension field, T is: ")
    T = EK.random_point()
    print(T)
    
    Z = Q + T
    print("Z is an addition of Q to the T")
    print(Z)


    # using the weil_pairing below
    # https://github.com/sagemath/sagesmc/blob/master/src/sage/schemes/elliptic_curves/ell_point.py#L1498
    x = Q.weil_pairing(T, 11178471414836345336059048855611390994257347541678490001885302105452889531220681887063760522435920242549498857070544319215681595199825520457978318663949088451956393179940626847131008230965896616721895736790011463601842300814881217357911666596763666142060364248293661581607975062108643224429446283393991255035656825117099140533934652875143386959036703603226566679927774482780478792984642946385754975582563741941170313319050311130987979716708470069385452052743295817756387189493033394289373456999033551211492693962602437835427803879042700605665271512616644853266461370443414567029631937208881495158644492620428589143183636655393718868986856362999803778403269108096669598360627088720895292022466546145805613907119407341002305902542929009712010952117502667466294030730061615881003093058836542115915779296546454823737486659498696382465941269092111594805775965304567737873324377973343180166011760630835377760362141121685804988129983376016160652759299142703468026467306556279618698252911728583946615756708459752059357920596901915230169254310895866402942405009235706175972318443692528574476593489954174142202606563337911501692785987709315325523637478127171293449780083027709715228710614005458956061504218741911709094828298418333340827389154912375715843217153929829606277277554688)
    y = T.weil_pairing(Q, 11178471414836345336059048855611390994257347541678490001885302105452889531220681887063760522435920242549498857070544319215681595199825520457978318663949088451956393179940626847131008230965896616721895736790011463601842300814881217357911666596763666142060364248293661581607975062108643224429446283393991255035656825117099140533934652875143386959036703603226566679927774482780478792984642946385754975582563741941170313319050311130987979716708470069385452052743295817756387189493033394289373456999033551211492693962602437835427803879042700605665271512616644853266461370443414567029631937208881495158644492620428589143183636655393718868986856362999803778403269108096669598360627088720895292022466546145805613907119407341002305902542929009712010952117502667466294030730061615881003093058836542115915779296546454823737486659498696382465941269092111594805775965304567737873324377973343180166011760630835377760362141121685804988129983376016160652759299142703468026467306556279618698252911728583946615756708459752059357920596901915230169254310895866402942405009235706175972318443692528574476593489954174142202606563337911501692785987709315325523637478127171293449780083027709715228710614005458956061504218741911709094828298418333340827389154912375715843217153929829606277277554688)
    z = Z.weil_pairing(Q, 11178471414836345336059048855611390994257347541678490001885302105452889531220681887063760522435920242549498857070544319215681595199825520457978318663949088451956393179940626847131008230965896616721895736790011463601842300814881217357911666596763666142060364248293661581607975062108643224429446283393991255035656825117099140533934652875143386959036703603226566679927774482780478792984642946385754975582563741941170313319050311130987979716708470069385452052743295817756387189493033394289373456999033551211492693962602437835427803879042700605665271512616644853266461370443414567029631937208881495158644492620428589143183636655393718868986856362999803778403269108096669598360627088720895292022466546145805613907119407341002305902542929009712010952117502667466294030730061615881003093058836542115915779296546454823737486659498696382465941269092111594805775965304567737873324377973343180166011760630835377760362141121685804988129983376016160652759299142703468026467306556279618698252911728583946615756708459752059357920596901915230169254310895866402942405009235706175972318443692528574476593489954174142202606563337911501692785987709315325523637478127171293449780083027709715228710614005458956061504218741911709094828298418333340827389154912375715843217153929829606277277554688)
    print("### The result of the weil pairing is ###")
    print(x)
    print(y)
    print(z)
    
    n = 11178471414836345336059048855611390994257347541678490001885302105452889531220681887063760522435920242549498857070544319215681595199825520457978318663949088451956393179940626847131008230965896616721895736790011463601842300814881217357911666596763666142060364248293661581607975062108643224429446283393991255035656825117099140533934652875143386959036703603226566679927774482780478792984642946385754975582563741941170313319050311130987979716708470069385452052743295817756387189493033394289373456999033551211492693962602437835427803879042700605665271512616644853266461370443414567029631937208881495158644492620428589143183636655393718868986856362999803778403269108096669598360627088720895292022466546145805613907119407341002305902542929009712010952117502667466294030730061615881003093058836542115915779296546454823737486659498696382465941269092111594805775965304567737873324377973343180166011760630835377760362141121685804988129983376016160652759299142703468026467306556279618698252911728583946615756708459752059357920596901915230169254310895866402942405009235706175972318443692528574476593489954174142202606563337911501692785987709315325523637478127171293449780083027709715228710614005458956061504218741911709094828298418333340827389154912375715843217153929829606277277554688
    
    U = int(order_P_Q/n**2)*EK.random_point()
    V = int(order_P_Q/n**2)*EK.random_point()
    print("THE ORDER OF U AND V ARE")
    #print(U.order())
    #print(V.order())


    # res is 1
    weil_u_v = U.weil_pairing(V, n)
    print(weil_u_v)



    #print("########## ABOUT Q ##############")
    #print(Q)
    #print("Q[0] is:    ")
    #print(Q[0])
    #print("Q[1] is:    ")
    #print(Q[1])
    #print("#### piP is")
    #piP = EK(Q[0]**q, Q[1]**q)
    #print(piP)
    #print(piP - Q)
    #orddd = 11178471414836345336059048855611390994257347541678490001885302105452889531220681887063760522435920242549498857070544319215681595199825520457978318663949088451956393179940626847131008230965896616721895736790011463601842300814881217357911666596763666142060364248293661581607975062108643224429446283393991255035656825117099140533934652875143386959036703603226566679927774482780478792984642946385754975582563741941170313319050311130987979716708470069385452052743295817756387189493033394289373456999033551211492693962602437835427803879042700605665271512616644853266461370443414567029631937208881495158644492620428589143183636655393718868986856362999803778403269108096669598360627088720895292022466546145805613907119407341002305902542929009712010952117502667466294030730061615881003093058836542115915779296546454823737486659498696382465941269092111594805775965304567737873324377973343180166011760630835377760362141121685804988129983376016160652759299142703468026467306556279618698252911728583946615756708459752059357920596901915230169254310895866402942405009235706175972318443692528574476593489954174142202606563337911501692785987709315325523637478127171293449780083027709715228710614005458956061504218741911709094828298418333340827389154912375715843217153929829606277277554688
    #Q.ate_pairing(T, orddd, 8, tr)

    
    E_D = EllipticCurve([(Fq4(a))/j**2, (Fq4(b))/j**3])
    print("D-twist\na=a/w^2={}\nb=b/w^3={}".format(E_D.a4(),E_D.a6()))
    print("PRINT FQ4")
    print(Fq4)
    N = [x for x in range(1000) if kronecker(x, q) == -1]; N
    print(N)
    print("The order of the twist is calculated...")
    order_twist = ZZ(105728290513165610986047842809751318841280324636548185919618761670766498489126319374080364643814474159560267765206481614470062820665803111802542081483890984079007411741658964681748397228590555002124666865154706881326247670437929905649005732882915063211059733358699029466871349150330576303159118955009483248411687688459258653412137668989908623154784117709434748451700616886945074033208009205147510709322838869285502070178802449805606834287959659872025482054893432653565073089295339866221304319849957430798002705359030400272401451744984691726540651157819242361547135100554273343028342409843532297688541312108594686251215425527313523033316)
    print(order_twist)
    assert order_twist % r == 0
    c2 = order_twist // r
    print("c2 is calculated..." )
    print(c2)
    print("check_curve_order(E_D, order_twist): {}".format(check_curve_order(E_D, order_twist)))
    
    exponent = (q**2 - q + 1)//r ## Is this correct for Lokum?
    e1 = exponent // q
    e0 = exponent % q
    
    e1_2naf = bits_2naf(e1)
    e0_2naf = bits_2naf(e0)
    print(e1_2naf)
    print(e0_2naf)
    
    ##

if __name__ == "__main__":
    print("Hello world!")

    arithmetic(False)
    #preparse("QQx.<x> = QQ[]")
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    # BLS12-377 seed
    # CHANGE: u0
    u0 = ZZ(0x8508C00000000001)
    # Lokum parameters
    r = ZZ(3618502788666131213697322783095070105623107215331596699973092056135872020481)
    q = ZZ(570227033427643938477351787526785557771164366858588273239167105706995178753794255646220989581954976296644199541099698255998850583874501725806067165976078609861)
    a = ZZ(408825162639581376425354251544709546275292709793540929634396921208238440519285149826351823244598728486772571705638177322608711025414928328755239670937630810045)
    b = ZZ(349808739615329235353531183724887241300117306968432494640082668760849877387756346906668495855253497480229722131408639028906859040081069711950364714123876057189)
    # cofactor of the curve
    c = ZZ(157586456811283429841422804972005719747057336671136400349968118675472488492656786156)
    tr= ZZ(3816706629471577761299121574721828830967739654627469433190085390763284409348826)
    D = ZZ(-143)

    Miller_loop_scalar = ZZ(3816706629471577761299121574721828830967739654627469433190085390763284409348825)
    assert Miller_loop_scalar == (tr-1)
    ## COMMENT FROM CP6 T = 506464946133393486072777102926336625944849939610982267859828541006717966526573193706126370441346337661774335955699621
    ## COMMENT FROM CP6 T = 0xcda91b419cd3f63a81dba4eb8fc81231d84f4491c1e1a4dfaf393e5c14beec4affbc93fba900003f8fbe3c000000007a5
    # absolute extension
    ## CHANGE? What's xi? I think it's the generator in Fq file? But, what's the generator?
    #xi = 13
    
    ## CHANGE? What's tr3
    #tr3=ZZ(-0x87219fb9ab6223dea08125f7b21a9d72bfd2eca93424662b01c6ba77d0efc9ada28e36672abc6ddb9ff90435350b0907a055422d2327d4b7662461c32b1a122e8cd4fdcbdcb496d5395acf3575580d5ef3b2859d2fb197730dcd15e9e306306baa548e380d51e48cf9f029557488a5d8e51a4d13319fa7b0ecf8d1613fd2fca09e981bf2d5f42d62a7ef0e68ce62bce5455ea)
    # we have r^2 divides q^6+1-t6 = (q3+1-t3)*(q^3+1+t3),
    # r divides (q+1-t) and (q+1-t) divides #E(Fq3)=(q3+1-t3)
    # r divides also the quadratic twist of E(Fq3) of order q^3+1+t3
    #order_twist = q^3 + 1 + tr3
    ## CHANGE? How to find the order_twist
    #order_twist = ZZ(0x2b87fda17134c96a73372daf9a1407c1bfe240a89db104c580df9d215037525cbc4ba24cf9ce08b93a5652a8e16a97fc618f319fb7fa85b5a9e180cd3f8dabebe38a09b261651d60f8f9f6fcd973093d8719c862c2bcd8dc3b1986ef18419cd0d1d6be027e0248d61e0c45d27076027ec67b23ca21268d5c2076b44f0a32c9a8abe9b50297df4c5efe04c0cec1019fb862c0494c868105fb739b5347526f005a82c4cb19e83b1a26fa18378c42661feec7f593c522696a885ca18780d64da536e647fec41f299502df523124f0b29870a703e5bde1bc49259646f3bbd388ab5642f80a80a89c72a4f78017ce6e9e8ccb996212d7771fa8a4fc5500b2352b94d58a035f2c4635f33ff5435f393055e21a030ae458054a1372030868e0e0540a3fca151d50b90)
    #assert (order_twist % r) == 0
    
    #c2=ZZ(0x19e70d3618ca3a1ad2ad61580ef78d7f8e3ef002559614ed25f7ebe6c44b3be7ce44dd084e05d4a86d8ebb7fdf52ad09033887d09b4d2c4b99b0d9b077d21a3b8a42c4e1e6f5a9eea59baa2aeb4a1fee44392814554748b589c61531dbce7d1d80da227dd39c1d8f9eda95d6abf3ef8c39cbf0df977800f25b044f4e9c9327ecaaf050acb91f3a597504275948226574956a099ace2c3c541018fe4d10cc7bedd37f629e8e6348eb48c92998f3358bb74eb9b6630846b7a5e400154f633cf4216f0595cee44d02cea7dc01aae2140556e918748f067931d6f388069c063c11f0c9a51593693f88c98a12bd486d2fb4b77fca151d50b90)
    #assert c2*r == order_twist
    # Final exponentiation
    # hard part: (q^2-q+1)/r = W0 + W1*q
    #W0 = ZZ(7000705447348627246181409558336018323010329260726930841638672011287206690002601216854775649561085256265269640040570922609783227469279331691880282815325569032149343779036142830666859805506518426649197067288711084398033)
    #W1 = ZZ(86482221941698704497288378992285180119495364068003923046442785886272123124361700722982503222189455144364945735564951562986)

    # easy final exp costs 4 + 34 + inversion + 18 = 56 + inversion ~= 81 m
    #W1_2naf=[0,1,0,1,0,-1,0,0,-1,0,-1,0,-1,0,1,0,0,0,1,0,1,0,0,-1,0,0,-1,0,0,0,0,-1,0,-1,0,-1,0,0,-1,0,0,1,0,-1,0,0,-1,0,0,0,0,-1,0,-1,0,1,0,-1,0,-1,0,-1,0,-1,0,0,0,0,-1,0,-1,0,-1,0,0,0,1,0,1,0,1,0,1,0,1,0,-1,0,0,0,-1,0,0,0,1,0,-1,0,0,0,0,-1,0,1,0,-1,0,1,0,1,0,-1,0,1,0,-1,0,0,0,0,1,0,1,0,0,0,0,-1,0,0,0,1,0,1,0,0,1,0,0,0,1,0,-1,0,1,0,-1,0,0,-1,0,-1,0,0,1,0,0,-1,0,0,1,0,1,0,1,0,0,0,0,1,0,0,0,0,1,0,0,1,0,1,0,-1,0,0,0,-1,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,-1,0,-1,0,0,0,1,0,1,0,1,0,0,0,0,0,-1,0,0,0,-1,0,1,0,-1,0,0,0,-1,0,1,0,-1,0,0,0,1,0,0,-1,0,0,0,0,-1,0,-1,0,-1,0,-1,0,1,0,1,0,0,1,0,0,0,0,-1,0,1,0,0,1,0,-1,0,-1,0,1,0,0,1,0,0,0,0,-1,0,0,-1,0,1,0,-1,0,-1,0,0,-1,0,1,0,-1,0,1,0,0,-1,0,0,0,0,0,0,-1,0,1,0,0,0,0,0,0,1,0,0,-1,0,1,0,0,0,-1,0,0,1,0,0,0,1,0,-1,0,-1,0,-1,0,0,0,0,1,0,1,0,0,0,-1,0,1,0,0,1,0,1,0,0,0,0,-1,0,0,0,0,0,1,0,0,1,0,-1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1,0,0,0,1]
    #W0_2naf=[1,0,0,0,1,0,-1,0,0,-1,0,-1,0,-1,0,1,0,0,0,1,0,1,0,0,-1,0,0,-1,0,0,0,0,-1,0,-1,0,-1,0,0,-1,0,0,1,0,-1,0,0,0,-1,0,0,0,-1,0,1,0,0,-1,0,-1,0,0,-1,0,-1,0,-1,0,0,0,0,0,1,0,0,0,1,0,0,1,0,1,0,0,1,0,-1,0,0,1,0,1,0,0,1,0,-1,0,0,0,-1,0,0,0,1,0,-1,0,0,-1,0,-1,0,0,-1,0,-1,0,0,-1,0,0,0,1,0,-1,0,-1,0,0,0,1,0,1,0,0,0,0,1,0,-1,0,0,1,0,0,1,0,1,0,0,1,0,0,0,0,0,-1,0,1,0,0,-1,0,0,0,1,0,-1,0,0,0,1,0,-1,0,0,0,0,0,0,0,0,-1,0,1,0,0,1,0,-1,0,0,1,0,1,0,0,1,0,0,0,0,1,0,1,0,0,-1,0,-1,0,0,1,0,0,0,0,-1,0,0,0,-1,0,-1,0,1,0,1,0,0,-1,0,1,0,-1,0,-1,0,0,0,0,1,0,1,0,0,-1,0,0,-1,0,0,-1,0,-1,0,0,-1,0,1,0,-1,0,-1,0,1,0,-1,0,0,1,0,0,0,0,1,0,-1,0,1,0,-1,0,-1,0,0,0,1,0,-1,0,0,0,0,1,0,0,0,1,0,0,0,-1,0,1,0,0,-1,0,0,-1,0,1,0,0,0,1,0,0,0,1,0,1,0,0,1,0,-1,0,0,0,-1,0,0,-1,0,1,0,0,1,0,1,0,1,0,0,0,-1,0,0,0,0,0,1,0,-1,0,0,1,0,0,0,1,0,1,0,0,0,1,0,-1,0,1,0,-1,0,0,0,1,0,0,0,1,0,-1,0,0,1,0,1,0,0,1,0,0,-1,0,0,0,0,-1,0,0,1,0,-1,0,-1,0,1,0,1,0,-1,0,0,0,-1,0,-1,0,0,-1,0,-1,0,1,0,0,0,1,0,-1,0,1,0,0,-1,0,1,0,-1,0,1,0,0,-1,0,0,1,0,1,0,0,-1,0,1,0,-1,0,1,0,0,0,0,0,0,1,0,-1,0,1,0,1,0,0,0,1,0,0,-1,0,1,0,-1,0,-1,0,0,-1,0,0,0,0,0,1,0,1,0,-1,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,1,0,-1,0,1,0,-1,0,1,0,1,0,1,0,-1,0,0,0,1,0,0,-1,0,1,0,0,-1,0,0,0,0,1,0,0,1,0,1,0,-1,0,0,-1,0,0,0,0,0,-1,0,1,0,-1,0,0,-1,0,0,0,-1,0,-1,0,0,1,0,-1,0,0,0,0,1,0,1,0,0,1,0,0,-1,0,-1,0,1,0,0,0,0,0,1,0,-1,0,-1,0,0,0,1,0,0,1,0,-1,0,0,0,-1,0,0,-1,0,-1,0,1,0,1,0,0,-1,0,0,0,0,-1,0,1,0,-1,0,0,0,0,0,0,0,0,0,1,0,-1,0,-1,0,1,0,-1,0,1,0,-1,0,0,1,0,0,0,0,1,0,0,0,0,-1,0,-1,0,-1,0,1,0,0,1,0,-1,0,-1,0,-1,0,-1,0,-1,0,0,1,0,0,1,0,0,1,0,-1,0,-1,0,0,0,1,0,1,0,0,0,1,0,1]

    #exponent = (q**2 - q + 1)//r
    #e1 = exponent // q
    #e0 = exponent % q
    #e1 += 1
    #e0 -= q
    #print("e1 = {}\ne0 = {}\nW1 = {}\nW0 = {}".format(e1, e0, W1, W0))
    #assert -W0 + W1*q == exponent
    #log2_e = float(log(abs(e1), float(2.0)) + log(abs(e0), float(2.0)))
    #log2_W = float(log(abs(W1), float(2.0)) + log(abs(W0), float(2.0)))
    #print("log_2(e1) + log_2(e0) = {}".format(log2_e))
    #print("log_2(W1) + log_2(W0) = {}".format(log2_W))
    #M = Matrix(ZZ, 2, 2, [exponent, 0, -q, 1])
    #R = M.LLL()
    #print(R)
    #for i in range(2):
    #    print("R[{}][0] + q*R[{}][1] = ({}) * ((q^2-q+1)//r) + ({})".format(i,i, (R[i][0] + R[i][1]*q) // exponent, (R[i][0] + R[i][1]*q) % exponent ))
    #log2_R0 = float(log(abs(R[0][0]), float(2.0)) + log(abs(R[0][1]), float(2.0)))
    #print("log_2(R00) + log_2(R01) = {}".format(log2_R0))
    #log2_R1 = float(log(abs(R[1][0]), float(2.0)) + log(abs(R[1][1]), float(2.0)))
    #print("log_2(R10) + log_2(R11) = {}".format(log2_R1))
    #assert W0 == -e0
    #assert W1 == e1
    test_lokum(a, b, q, r, tr, c, D)