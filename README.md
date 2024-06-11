# projekt1_poprawione
# Instrukcja obsługi programu skrypt.py

Program skrypt.py umożliwia przekształcanie współrzędnych geodezyjnych między różnymi układami odniesienia oraz formatami, korzystając z dostępnych elipsoid geodezyjnych.

## Wymagania:

Program jest przeznaczony dla systemu operacyjnego Windows i wymaga zainstalowanego Pythona w wersji 3.9 oraz następujących bibliotek:
- numpy (importowane jako np): używane do operacji na macierzach i wektorach.
- math: zawiera podstawowe funkcje matematyczne.
- sys: zawiera metody i zmienne służące do modyfikowania wielu elementów środowiska wykonawczego języka Python.

## Elipsoidy geodezyjne:

Program obsługuje następujące elipsoidy geodezyjne:
- WGS84: Elipsoida zdefiniowana w ramach systemu WGS84.
- GRS80: Elipsoida zdefiniowana w ramach systemu GRS80.
- mars: Elipsoida geodezyjna Marsa.

## Dostępne transformacje:

1. XYZ_BLH: Przekształcenie współrzędnych XYZ na współrzędne geodezyjne BLH.
2. BLH_XYZ: Przekształcenie współrzędnych BLH na współrzędne XYZ.
3. XYZ_NEU: Przekształcenie współrzędnych XYZ na współrzędne NEU.
4. BL_PL2000: Przekształcenie współrzędnych BLH na współrzędne PL-2000.
5. BL_PL1992: Przekształcenie współrzędnych BLH na współrzędne PL-1992.

## Sposób użycia:

Program uruchamia się poprzez podanie parametrów przez wiersz poleceń.

## Przykłady użycia:

1. Przekształcenie współrzędnych XYZ na współrzędne geodezyjne BLH:
```
python skrypt.py -plik dane.txt -elip WGS84 -funkcja XYZ_BLH
```

2. Przekształcenie współrzędnych BLH na współrzędne XYZ:
```
python skrypt.py -plik dane.txt -elip WGS84 -funkcja BLH_XYZ
```

3. Przekształcenie współrzędnych XYZ na współrzędne NEU:
```
python skrypt.py -plik dane.txt -elip WGS84 -funkcja XYZ_NEU
```

4. Przekształcenie współrzędnych BLH na współrzędne PL-2000:
```
python skrypt.py -plik dane.txt -elip WGS84 -funkcja BL_PL2000
```

5. Przekształcenie współrzędnych BLH na współrzędne PL-1992:
```
python skrypt.py -plik dane.txt -elip WGS84 -funkcja BL_PL1992
```

## Parametry:

- -plik: Nazwa pliku z danymi wejściowymi (pełna nazwa z rozszerzeniem, np. "dane.txt").
- -elip: Wybór elipsoidy (dostępne wartości: "WGS84", "GRS80", "mars").
- -funkcja: Wybór rodzaju transformacji (dostępne wartości: "XYZ_BLH", "BLH_XYZ", "XYZ_NEU", "BL_PL2000", "BL_PL1992").

## Format danych wejściowych:

Dane wejściowe powinny być przechowywane w pliku tekstowym w tym samym folderze co program skrypt.py. Współrzędne XYZ lub BLH powinny być oddzielone przecinkami. Jednostka danych wejściowych zależy od wybranej funkcji.

## Format danych wyjściowych:

Wyniki przekształcenia zostaną zapisane w pliku tekstowym o nazwie "WYNIK_NAZWA_FUNKCJI.txt", gdzie "NAZWA_FUNKCJI" to nazwa funkcji transformacji (np. "WYNIK_XYZ_BLH.txt"). Współrzędne będą oddzielone spacjami, a każda linia będzie zawierać współrzędne dla jednego punktu. Jednostka danych wyjściowych zależy od wybranej funkcji.

## Błędy:

Program obsługuje następujące błędy:
- Błąd podania niewłaściwej liczby argumentów.
- Błąd przy podaniu parametrów bez wartości.
- Błąd przy podaniu nieprawidłowego modelu elipsoidy.
- Błąd przy podaniu nieobsługiwanej transformacji.
- Błąd w przypadku braku pliku wejściowego.
- Błąd w przypadku pliku o niewłaściwym formacie.

Program rozpoczyna wczytywanie danych po pierwszych czterech linijkach nagłówka pliku wejściowego.
