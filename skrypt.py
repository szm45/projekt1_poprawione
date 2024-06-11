
import numpy as np
from argparse import ArgumentParser
import sys

class Transformacje:
    
    
    def __init__(self, elipsoida):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: https://en.wikibooks.org/wiki/PROJ.4#Spheroid
        + Parametry planet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
        """
        self.a = elipsoida[0]
        self.e2 = elipsoida[1]
    def dms(self, x):
        '''
        Funkcja dms służy nam do zamiany jednostek, z radianów na stopnie.  

        Parametry
        ----------
        x : FLOAT
            [rad].

        Returns
        x : STR
            [dms] - stopnie, minuty, sekundy
        '''
        sig = ' '
        if x < 0:
            sig = '-'
            x = abs(x)
        x = x * 180/np.pi
        d = int(x)
        m = int(60 * (x - d))
        s = (x - d - m/60)*3600
        if s > 59.999995:
            s = 0
            m = m + 1
        if m == 60:
            m = 0
            d = d + 1
        
        d = str(d)
        if len(d) == 1:
            d = "  " + d
        elif len(d) == 2:
            d = " " + d
        elif len(d) == 3:
            d = d
            
        if m < 10:
            m = str(m)
            m = "0" + m
            
        if s < 10:
            s = "%.5f"%s
            s = str(s)
            s= "0" + s
        else:
            s = ("%.5f"%s)
            
        x1=(f'{d}°{m}′{s}″')  
        return(x1)
        
        
    def NP(self, fi):
        '''
       Funkcja oblicza promień przekroju w pierwszej pionowej płaszczyźnie, który jest potrzebny m.in. do zastosowania algorytmu Hirvonena
        
        Parameters
        ----------
        fi : FLOAT
            [radiany] - szerokość geodezyjna

        Returns
        -------
        N : float
            [metry] - promień przekroju w pierwszym wertykale

        '''
        N = self.a / np.sqrt(1 - self.e2 * np.sin(fi)**2)
        return(N)
    
    
    def algorytm_hirvonena(self, X, Y, Z, output="dec_degree"):
        '''
        Algorytm Hirvonena służy do konwersji współrzędnych prostokątnych (x, y, z) na współrzędne geodezyjne (fi, Lambda, h). Jest to proces iteracyjny, który umożliwia osiągnięcie dokładności rzędu 1 milimetra po kilku powtórzeniach procedury.

         Parametry
         ----------
         X, Y, Z : FLOAT
              współrzędne w układzie orto-kartezjańskim, 

         Returns
         -------
         fi : FLOAT
         [st. dz.] - szerokość geodezyjna.
         lam : FLOAT
         [st. dz.] - długośc geodezyjna.
         h : FLOAT
         [m] - wysokość elipsoidalna
         output [STR] - opcjonalne, domylne 
         -dec_degree - st. dz
         -dms - st, min, sek
         -radiany - rad 
         """
        '''

        p = np.sqrt(X**2 + Y**2)
        fi = np.arctan(Z/(p * (1 - self.e2)))
        while True:
            N = Transformacje.NP(self, fi)
            h = (p / np.cos(fi)) - N
            fip = fi
            fi = np.arctan(Z / (p * (1 - self.e2 * (N / (N+h)))))
            if np.abs(fip - fi) < (0.000001/206265):
                break
        lam = np.arctan2(Y, X)
        if output == "dec_degree":
            fi=(fi*180/np.pi)
            lam=(lam*180/np.pi)
            return (fi, lam, h)
        elif output == "dms":
            fi = Transformacje.dms(self, fi)
            lam = Transformacje.dms(self, lam)
            return (fi,lam,h) 
        elif output == 'radiany':
            
            
            return(fi,lam,h)
        else:
            raise NotImplementedError(f"{output} - output format not defined")

    def odwrotny_hirvonen(self, fi, lam, h):
       '''
      Algorytm odwrotny do algorytmu Hirvonena służy do przekształcania współrzędnych geodezyjnych (B, L, H) 
      na współrzędne ortokartezjańskie (x, y, z). To proces umożliwiający przejście z opisu geodezyjnego punktu 
      na powierzchni ziemi do odpowiadającej mu lokalizacji w trójwymiarowym układzie współrzędnych kartezjańskich.
       Parametry
       ----------
       fi : FLOAT
           [st. dz] - szerokość
       lam : FLOAT
           [st. dz] - długośc 
       h : FLOAT
           [m] - wysokość 
       Returns
       -------
        X, Y, Z : FLOAT
             [m] - współrzędne orto-kartezjańskie

       '''
       fi=np.radians(fi)
       lam=np.radians(lam)
       
       N = Transformacje.NP(self, fi)
       X = (N + h) * np.cos(fi) * np.cos(lam)
       Y = (N + h) * np.cos(fi) * np.sin(lam)
       Z = (N *(1-self.e2) + h) * np.sin(fi) 
       
       return(X,Y,Z)
   
   
    def flh2PL1992(self, fi, lam):
       '''
       Układ współrzędnych 1992 (PUWG-92) to system płaskich współrzędnych prostokątnych,
       który używa odwzorowania Gaussa-Krügera dla elipsoidy GRS80 w ramach pojedynczej dziesięciostopniowej strefy.

       Parametry
       ----------
       fi : FLOAT
           [st. dz.] - szerokość
       lam : FLOAT
           [st. dz] - długośc

       Returns
       -------
        X1992, Y1992 : FLOAT
             [m] - współrzędne (1992)

       '''
       
       if lam > 25.5 or lam < 13.5:
           raise NotImplementedError(f"{Transformacje.dms(self, np.radians(lam))} ten południk nie jest obsługiwany przez układ współrzędnych płaskich PL1992")
           
       if fi > 55 or fi < 48.9:
           raise NotImplementedError(f"{Transformacje.dms(self, np.radians(fi))} ten równoleżnik nie jest obsługiwany przez układ współrzędnych płaskich PL1992")
           
       fi = np.radians(fi)
       lam = np.radians(lam)
       a2 = self.a**2
       b2 = a2 * (1 - self.e2)
       e_2 = (a2 - b2)/b2
       l0 = np.radians(19)
       dl = lam - l0
       dl2 = dl**2
       dl4 = dl**4
       t = np.tan(fi)
       t2 = t**2
       t4 = t**4
       n2 = e_2 * (np.cos(fi)**2)
       n4 = n2 ** 2
       N = Transformacje.NP(self, fi)
       e4 = self.e2**2
       e6 = self.e2**3
       A0 = 1 - (self.e2/4) - ((3*e4)/64) - ((5*e6)/256)
       A2 = (3/8) * (self.e2 + e4/4 + (15*e6)/128)
       A4 = (15/256) * (e4 + (3*e6)/4)
       A6 = (35*e6)/3072
       sigma = self.a * ((A0 * fi) - A2 * np.sin(2*fi) + A4 * np.sin(4*fi) - A6 * np.sin(6*fi))
       xgk = sigma + ((dl**2)/2) * N * np.sin(fi) * np.cos(fi) * (1 + ((dl**2)/12)*(np.cos(fi)**2)*(5 - t2 + 9 * n2 + 4 * n4) + (dl4/360) * (np.cos(fi)**4)*(61 - (58 * t2) + t4 + (270 * n2) - (330 * n2 * t2)))
       ygk = dl * N * np.cos(fi) * (1 + (dl2/6) * (np.cos(fi)**2) * (1 - t2 + n2) + (dl4/120) * (np.cos(fi)**4) * (5 - (18 * t2) + t4 + (14 * n2) - 58 * n2 * t2))
       x92 = xgk * 0.9993 - 5300000
       y92 = ygk * 0.9993 + 500000
       return(x92,y92)
   
   
    def flh2PL2000(self, fi, lam):
       '''
      Układ współrzędnych 2000 to system prostych współrzędnych płaskich.
      Wykorzystuje on odwzorowanie Gaussa-Krügera dla elipsoidy GRS 80 w czterech określonych strefach, na południkach
      15°E, 18°E, 21°E i 24°E.

       Parametry
       ----------
       fi : FLOAT
           [st. dz.] - szerokość 
       lam : FLOAT
           [st. dz.] - długośc 

       Returns
       -------
        X2000, Y2000 : FLOAT
             [m] - współrzędne 

       '''
         
       if lam >= 13.5 and lam < 16.5:
           l0 = np.radians(15)
       elif lam >= 16.5 and lam < 19.5:
           l0 = np.radians(18)
       elif lam >= 19.5 and lam < 22.5:
           l0 = np.radians(21)
       elif lam >= 22.5 and lam <= 25.5:
           l0 = np.radians(24)
       else:
           raise NotImplementedError(f"{Transformacje.dms(self, np.radians(lam))} ten południk nie mieści się w zakresie")
       
       if fi > 55 or fi < 48.9:
           raise NotImplementedError(f"{Transformacje.dms(self, np.radians(fi))} ten równoleżnik nie mieści się w zakresie")
           
       fi = np.radians(fi)
       lam = np.radians(lam)
       a2 = self.a**2
       b2 = a2 * (1 - self.e2)
       e_2 = (a2 - b2)/b2
       dl = lam - l0
       dl2 = dl**2
       dl4 = dl**4
       t = np.tan(fi)
       t2 = t**2
       t4 = t**4
       n2 = e_2 * (np.cos(fi)**2)
       n4 = n2 ** 2
       N = Transformacje.NP(self, fi)
       e4 = self.e2**2
       e6 = self.e2**3
       A0 = 1 - (self.e2/4) - ((3*e4)/64) - ((5*e6)/256)
       A2 = (3/8) * (self.e2 + e4/4 + (15*e6)/128)
       A4 = (15/256) * (e4 + (3*e6)/4)
       A6 = (35*e6)/3072
       sigma = self.a * ((A0 * fi) - A2 * np.sin(2*fi) + A4 * np.sin(4*fi) - A6 * np.sin(6*fi))
       xgk = sigma + ((dl**2)/2) * N * np.sin(fi) * np.cos(fi) * (1 + ((dl**2)/12)*(np.cos(fi)**2)*(5 - t2 + 9 * n2 + 4 * n4) + (dl4/360) * (np.cos(fi)**4)*(61 - (58 * t2) + t4 + (270 * n2) - (330 * n2 * t2)))
       ygk = dl * N * np.cos(fi) * (1 + (dl2/6) * (np.cos(fi)**2) * (1 - t2 + n2) + (dl4/120) * (np.cos(fi)**4) * (5 - (18 * t2) + t4 + (14 * n2) - 58 * n2 * t2))
       strefa = round(l0 * 180/np.pi)/3
       x00 = xgk * 0.999923
       y00 = ygk * 0.999923 + strefa * 1000000 + 500000
       return(x00,y00)
   
    def dXYZ(self, xa, ya, za, xb, yb, zb):
        '''
       Funkcja ta służy do oblicenia różnic, pomiędzy wspolrzednymi, punktów A oraz B.

        Parametry
        ----------
        XA, YA, ZA, XB, YB, ZB: FLOAT
             wsp. orto-kart. 

        Returns
        -------
        dXYZ : ARRAY
            macierz różnicy współrzędnych

        '''
        dXYZ = np.array([xb-xa, yb-ya, zb-za])
        return(dXYZ)
    
    
    def rneu(self, fi, lam):
        '''
        Funkcja ta służy do tworznia macierzy obrotu (R),
        Macierz ta jest konieczna do macierzy NEU

        Parametry
        ----------
        fi : FLOAT
            [st. dz.] - szerokość 
        lam : FLOAT
            [st. dz.] - długośc 

        Returns
        -------
        R ARRAY
            macierz obrotu R
             
        '''
        fi=np.radians(fi)
        lam=np.radians(lam)
        R = np.array([[-np.sin(fi)*np.cos(lam), -np.sin(lam), np.cos(fi)*np.cos(lam)],
                      [-np.sin(fi)*np.sin(lam),  np.cos(lam), np.cos(fi)*np.sin(lam)],
                      [np.cos(fi),             0,         np.sin(fi)          ]])
        return(R)
    
    
    def xyz2neu(self, fi, lam, xa, ya, za, xb, yb, zb):
        '''
        Układ współrzędnych horyzontalnych opisuje położenie obiektów astronomicznych względem
        lokalnego horyzontu. Zenit i nadir, są ważnymi punktami w tym układzie. Współrzędne horyzontalne
        zmieniają się wraz z ruchem obserwatora i czasem, co pozwala określić chwilową pozycję obiektów na niebie.

        Parametry
        ----------
        fi : FLOAT
            [st. dz.] - szerokość 
        lam : FLOAT
            [st. dz.] - długośc
        XA, YA, ZA, XB, YB, ZB: FLOAT
             współrzędne orto-kartezjańskie

        Returns
        -------
        neu : STR
            wsp. horyz.
            

        '''
        dX = Transformacje.dXYZ(self, xa, ya, za, xb, yb, zb)
        R = Transformacje.rneu(self, fi,lam)
        neu = R.T @ dX
        n = neu[0];   e = neu[1];   u = neu[2]
        n = "%.16f"%n; e = "%.16f"%e; u="%.16f"%u
        dlugosc = []
        xx = len(n); dlugosc.append(xx)
        yy = len(e); dlugosc.append(yy)
        zz = len(u); dlugosc.append(zz)
        P = 50
        
        while xx < P:
            n = str(" ") + n
            xx += 1
        
        while yy < P:
            e = str(" ") + e
            yy += 1
            
        while zz < P:
            u = str(" ") + u
            zz +=1
            
        return(n, e, u)
    
    def wczytanie(self, Dane):
       with open (Dane,"r") as plik:
           tip = np.genfromtxt(plik, delimiter=",",dtype = '<U20', skip_header = 4)
           X=[]
           Y=[]
           Z=[]
           for i in tip:
               x=i[0]
               X.append(float(x))
               y=i[1]
               Y.append(float(y))
               z=i[2]
               Z.append(float(z))
           wielkosc = len(X)
       return(X, Y, Z, wielkosc)
    
    def zapisanie(self, X, Y, Z, f, l, h, x92, y92, x00, y00, N, E, U, xyz_txt, neu_txt ): 
        '''
        funkcja ta służy do zapania wynikiów obliczeń (x, y, z, f, l, h, x92, y92, x1992, y1992, x2000, y2000 ,neu).
        Następnie, z uskanych wynów tworzy tabele.
    
        Parametry
        ----------
        XYZ : LIST
             [m] - współrzędne w układzie orto-kartezjańskim, 
         fi : LIST
             [dms] - szerokość
         lam : LIST
             [dms] - długośc 
         h : LIST
             [m] - wysokość elipsoidalna
        X1992, Y1992 : LIST
             [m] - współrzędne w PL1992
         X2000, Y2000 : LIST
             [m] - współrzędne w PL2000
        neu : str
            wsp. horyz.
    
        Returns
        -------
        PLIK TXT
    
        '''
        for i in range(len(X)):
            X[i] = Transformacje.zmiana_na_dms(self, X[i])
            Y[i] = Transformacje.zmiana_na_dms(self, Y[i])
            Z[i] = Transformacje.zmiana_na_dms(self, Z[i])
        
        with open(xyz_txt , "w",  encoding="utf-8") as plik:
            plik.write(f"Wyniki_obliczen; XYZ, fi, lambda, h, x1992, y1992, x2000, y2000.\n")
            plik.write(f"Znak '-' w; x1992, y1992, x2000, y2000 oznacza, że dla podanych współrzędnych (X, Y, Z) po obliczeniu współrzędnych geodezyjnych fi oraz lam. Fi i lam nie należą do dozwolonych współrzędnych układów 1992 oraz 2000.\n")
            plik.write(f"\n")
            plik.write(f"          X                    Y                    Z                    fi                 lambda                 h                  x1992                y1992                x2000                y2000        ")
            plik.write(f"\n")
            for x, y, z, f, l, h, x92, y92, x00, y00 in zip(X, Y, Z, f, l, h, x92, y92, x00, y00):
                plik.write(f"{x}{y}{z}     {f}     {l}{h}{x92}{y92}{x00}{y00}")
                plik.write(f"\n")
        
        with open(neu_txt , "w", encoding="utf-8") as plik1:
            plik1.write(f"Wyniki_obliczen; neu.\n")
            plik1.write(f"                        n                                                 e                                                 u                         ")
            plik1.write(f"\n")
            for n, e, u in zip(N, E, U):
                plik1.write(f"{n}{e}{u}")
                plik1.write(f"\n")
        
    def wczytanie_oraz_zapisanie(self, Dane, output ='dms' , xyz_txt = 'wszystkie_wyniki.txt', neu_txt = "Wyniki_NEU.txt" ):
        '''
         funkcja ta wczytuje i zapisuje plik.
    
         Parameters
         ----------
         Dane : txt
             Plik z danymi xyz.
         output : str
             sposób w jakiej ma zapisywać współrzędne fi, lam [dms, rad, st. dzies.] .
         XYZ_txt: STR
             nazwa pliku wynikowego na xyz, flh, PL1992, PL2000
         NEU_txt: STR
    
         Returns
         -------
         Plik txt
    
         '''
        X, Y, Z, C = Transformacje.wczytanie(self, Dane)
        F=[]
        L=[]
        H=[]
        X92=[]
        Y92=[]
        X00=[]
        Y00=[]
        N=[]
        E=[]
        U=[]
        for x, y, z in zip(X, Y, Z):
             f,l,h = Transformacje.algorytm_hirvonena(self, x, y, z, output = output)
             if output == "dms":
                 F.append(f)
                 L.append(l)
             elif output == "radiany":
                 f=Transformacje.zmiana_na_rad(self,f)
                 l=Transformacje.zmiana_na_rad(self,l)
                 F.append(f)
                 L.append(l)
             else:
                 f=Transformacje.zmiana_na_dziesietne(self,f)
                 l=Transformacje.zmiana_na_dziesietne(self,l)
                 F.append(f)
                 L.append(l)
             H.append(Transformacje.zmiana_na_dms(self, h))
             f,l,h = Transformacje.algorytm_hirvonena(self, x, y, z)
             
             if l >= 13.5 and l <= 25.5 and f <= 55.0 and f >= 48.9:
                 x92, y92 = Transformacje.flh2PL1992(self, f,l)
                 X92.append(Transformacje.zmiana_na_dms(self, x92))
                 Y92.append(Transformacje.zmiana_na_dms(self, y92))
                 x00, y00 = Transformacje.flh2PL2000(self, f,l)
                 X00.append(Transformacje.zmiana_na_dms(self, x00))
                 Y00.append(Transformacje.zmiana_na_dms(self, y00))
             else:
                 x92 = "         '-'         " ; X92.append(x92)
                 y92 = "         '-'         " ; Y92.append(y92)
                 x00 = "         '-'         " ; X00.append(x00)
                 y00 = "         '-'         " ; Y00.append(y00)
         
        f1, l1, h1 = Transformacje.algorytm_hirvonena(self, X[0], Y[0], Z[0])
        n1, e1, u1 = Transformacje.xyz2neu(self, f1, l1, X[0], Y[0], Z[0], X[-1], Y[-1], Z[-1])
        N.append(n1)
        E.append(e1)
        U.append(u1)
         
        i=0
        while i<(C-1):
            f, l, h = Transformacje.algorytm_hirvonena(self, X[i], Y[i], Z[i])
            n, e, u = Transformacje.xyz2neu(self, f, l, X[i], Y[i], Z[i], X[i+1], Y[i+1], Z[i+1])
            N.append(n)
            E.append(e)
            U.append(u)
            i+=1
         
        Transformacje.zapisanie(self, X, Y, Z, F, L, H, X92, Y92, X00, Y00, N, E, U, xyz_txt, neu_txt )     
         
         
    def zmiana_na_rad(self, liczba):
        '''
        zamienia z float na str, zachowuje odpowiednią dokładnosc jednostek przy wynikach końcowych.

             Parameters
        ----------
        liczba : FLOAT
            Wynik podany w rad.
             Returns
        -------
        liczba : STR
            string
             '''
        zm_liczba = "%.12f"%liczba
        P = 16
        xx = len(zm_liczba)
        while xx < P:
            zm_liczba = str(" ") + zm_liczba
            xx += 1
        return(zm_liczba)  
         
    def zmiana_na_dziesietne(self, liczba):
        '''
        zamienia z float na str z odpowiednią dokładnocia jednostek dzies, przy wyniku koncowym.
        ----------
        liczba : FLOAT
            Wynik podany w st. dzies.
             Returns
        -------
        liczba : STR
            string
             '''
        zm_liczba = "%.10f"%liczba
        P = 16
        xx = len(zm_liczba)
        while xx < P:
            zm_liczba = str(" ") + zm_liczba
            xx += 1
        return(zm_liczba)        
        
    
    def zmiana_na_dms(self, liczba):
        '''
        funkcja ta wpływa na zapis wyniku końcowego. 
             Parameters
        ----------
        liczba : FLOAT
            Wynik podane w DMS.
             Returns
        -------
        liczbe : STR
            string
        
        '''
        zm_liczba = "%.3f"%liczba
        P = 21
        xx = len(zm_liczba)
        while xx < P:
            zm_liczba = str(" ") + zm_liczba
            xx += 1
        return(zm_liczba)
           
if __name__ == "__main__":
    elipsoidy = {
        'WGS84': [6378137.000, 0.00669437999013],
        'GRS80': [6378137.000, 0.00669438002290],
        'Krasowski': [6378245.000, 0.00669342162296]
    }

    funkcje = {
        'XYZ_BLH': 'xyz2flh',
        'BLH_XYZ': 'flh2xyz',
        'XYZ_NEU': 'xyz2neu',
        'BL_PL2000': 'PL2000',
        'BL_PL1992': 'PL1992'
    }

    argumenty = sys.argv[1:]

    if not all(param in argumenty for param in ['-plik', '-elip', '-funkcja']):
        raise Exception('Nie podano wszystkich argumentów: -plik, -elip, -funkcja. Podaj wszystkie argumenty.')

    try:
        elip = argumenty[argumenty.index('-elip') + 1]
        trans_wsp = argumenty[argumenty.index('-funkcja') + 1]
        plik = argumenty[argumenty.index('-plik') + 1]
    except IndexError:
        raise Exception('Nie podano wartości dla wszystkich wymaganych parametrów')

    if elip not in elipsoidy:
        raise Exception('Podano niewłaściwy typ modelu elipsoidy. Podaj jeden z dostępnych: GRS80, WGS84, Krasowski.')

    elipsoida = elipsoidy[elip]
    geo = Transformacje(elipsoida)

    if trans_wsp not in funkcje:
        raise Exception('Skrypt nie obsługuje podanej transformacji. Podaj jedną z możliwych: XYZ_BLH, BLH_XYZ, XYZ_NEU, BL_PL2000, BL_PL1992.')

    try:
        with open(plik, 'r') as f, open(f"WYNIK_{trans_wsp.upper()}.txt", 'w') as wynik:
            linie = f.readlines()[4:]  # Pominięcie pierwszych czterech linii
            for index, linia in enumerate(linie):
                linia = linia.strip()
                if trans_wsp in ['XYZ_BLH', 'BLH_XYZ', 'XYZ_NEU']:
                    x_str, y_str, z_str = linia.split(',')
                    x, y, z = float(x_str), float(y_str), float(z_str)
                    if trans_wsp == 'XYZ_BLH':
                        result = geo.algorytm_hirvonena(x, y, z)
                    elif trans_wsp == 'BLH_XYZ':
                        result = geo.odwrotny_hirvonen(x, y, z)
                    elif trans_wsp == 'XYZ_NEU':
                        if index == 0:
                            X0, Y0, Z0 = x, y, z
                            continue
                        result = geo.xyz2neu(X0, Y0, Z0, x, y, z)
                    wynik.write(' '.join(map(str, result)) + '\n')
                else:
                    fi_str, lam_str = linia.split(',')
                    fi, lam = float(fi_str), float(lam_str)
                    if trans_wsp == 'BL_PL2000':
                        result = geo.flh2PL2000(fi, lam)
                    elif trans_wsp == 'BL_PL1992':
                        result = geo.flh2PL1992(fi, lam)
                    wynik.write(' '.join(map(str, result)) + '\n')
    except FileNotFoundError:
        raise Exception('Podany plik nie istnieje. Podaj inny plik, sprawdź jego lokalizację lub sprawdź nazwę podanego pliku.')
    except (KeyError, IndexError, ValueError):
        raise Exception('Nieodpowiedni format pliku.')
    except AttributeError:
        print("Podana funkcja/elipsoida nie istnieje, proszę wprowadzić dostępne wartości.")

    print('Zapisano. Wyniki znajdują się w pliku WYNIK_<funkcja>.txt')