NAMI-projektin R-skriptien versio 0.5.8
Jarno Tuimala, 3.10.2008

Sisällys:
- Tehtävää
- Versioon kuuluvat skriptit
- Nykyiset ominaisuudet
- Muutokset edellisestä versiosta
- R:n asennus Coronalle (tarvittavat paketit)

---------------------------------------------------------------------------------------------------------------------------

Tehtävää:
- Bugeja:
   - yeast2-siruilla ei toimi hypergeometrinen testi, koska ENTREZID:tä ei löydy annotaatiokirjastoista.
      - use systematic names
   - väritä PCA/NMDS-kuvat group-sarakkeen perusteella
   - volcanoplot skriptin pitäisi palauttaa myös data
   - Tee annotaatiokirjastot:
      - Illumina mouse V2
      - Illumina human V3
      - Agilentin miRNA-sirut
   - norm-cdna ja norm-agilent eivät anna flagejä, koska niiden rivimäärä eroaa datan rivimäärästä
     (koska datasta on laskettu replikaattispottien suhteen keskiarvo, mutta flageistä ei)
   - Kuvien antamat tulokset myös tekstinä
	- Kaikki tool-sivut ajan tasalle
         - Supported-chiptype -sivulle myös sirujen nimet systeemissä (R:ssä esim hgu133a)
         - artikkeliviitteet
	- Getting started -sivu
- Kehitettävää:
   - skriptit
      - Fold changen klusterointi
         - klusterointi-skripteihin valinta [chip., FC.] tms, jossa voi valita
           haluaako klusteroida chip. vai FC. sarakkeet.
      - Pathway analysis: show which genes contributes to the results
      - Ability to tell how big GO classes are taken into account in the analysis
      - Gene set test: show y-axis values
      - Arabidopsis promoter analysis
      - Many probes for the same gene: combine the expression values to one entry
      - kuvaavammat parametrinimet linear modellingissa (esim. main.effect1, technical.replicate, pairing, main.effect1.is.factor) 
      - storeyn Q-arvokorjaus yhdeksi FDR-vaihtoehdoksi
      - sort data 
      - extract rows
      - merge a table and a genelist
      - laajenna search-by-genename toiminnallisuutta
      - Bring in own annotations
         - New tools for running analyses using these
      - general normalization
          - mahdollisuus määritellä pudotusvalikoilla käytettävät sarakkeet
          - Agilentille oma skripti siten, että tab-delimited fileet vain tuodaan sisään,
            ja skripti hoitaa kaiken muun?
      - HGU95A-sarjan siruja custom.chiptype:ksi norm-affy:yyn
      - Rikastuneiden GO-luokkien analyysin pitäisi palauttaa geenitkin (miten?)
      - MA-plot for Illumina
      - direct data import form GeneSpring
	- survival analysis -skripti
      - extract genes using a category (GO, KEGG) term
      - search by gene ontology class
      - Normalisointi positiivisten kontrollien (spiked controls) perusteella
      - normalisointi tiettyihin geeneihin
      - plot residuals from RLE/NUSE
      - permutation for multiple testing correction
      - Help in selecting a correct p-value cutoff (Qvalue package)
      - connectivity map?
      - other ontologies (cMAP, BioCarta, reactome, panther)
	   - using BioPAX through Rredland
      - Fix export to soft to include all information from phenodata
   - more wizards
      - GO-classification based on min 3-fold change
      - extract the significant pathways
      - QC
   - ISOT JUTUT
      - Gene arrays
      - Tiling arrays
      - ChIP-on-chip
      - CGH
      - SNP analysis
         - add support for Illumina chips
         - link to Plink
         - Sampsa Hautaniemen ryhmän tuotokset
            - CohortComparator
               - Tsekkaa fortran-koodin kääntäminen R:ää varten
	- fylogenetiikka
	- geenikartoitus
- R-kirjasto-ongelmia:
   - Seuraavat kirjastot eivät vielä asentuneet:
	- Rgraphviz
   - Skriptit, jotka eivät siksi toimi:
	- pathway-bayesian.R

---------------------------------------------------------------------------------------------------------------------------

Tähän versioon kuuluvat skriptit:



---------------------------------------------------------------------------------------------------------------------------

Nykyiset ominaisuudet:



---------------------------------------------------------------------------------------------------------------------------

Muutokset edellisestä versiosta (0.5.7)
- Calculate-fold-change toimii nyt myös vain kahdella sarakkeella. Aiemmin skripti kaatui siihen, että kun
  kahden sarakkeen data framesta poimitaan yksi sarake, se konvertoidaan automaattisesti vektoriksi, jolle
  ncol()-komento antaa dimensioiksi NULL. Homma korjattu siten, että vektori muutetaan aina ensin data frameksi
  ennen laskentaa.
- Hierarkkisen klusteroinnin palauttama tulostiedosto muutettu hc.txt:stä hc.tre:ksi yhteensopivuusongelmien
  välttämiseksi.
- change-interpretation -skriptissä ollut erhe, joka jätti flag-sarakkeet dataan (ja kaatoi siis skriptin) on 
  korjattu.
- plot-dendrogram-skriptin oletusarvoksi muutettu chips genes sijaan.
- Uusi skripti: Laura Elo-Uhlgrenin ROTS
- Modifioitu normalisointi-tyokaluja siten, että normalisoituun dataan kirjoitetaan SYMBOL- ja DESCRIPTION-sarakkeet.
- calculate-fold-change -skripti antaa nyt myös datan FC:n lisäksi. FD muutettu tulosteessa muotoon FC.
- stat-two-groups: FD muutettu muotoon FC.
- stat-timeseries ei näyttänyt parametria p.value.threshold. Tämä johtui siitä, että edellisen parametrin perästä
  puuttuivat sulut. Sulut lisätty.
- norm-agilent- ja norm-cdna-skripteihin lisätty between arrays -normalisointioptioksi Aquantile.
- norm-affy -skriptin pitäisi nyt mennä läpi, vaikkei annotaatiokirjastoa löytyisi, jolloin SYMBOL ja DESCRIPTION-
  sarakkeita ei dataan kirjoiteta. 


Muutokset edellisestä versiosta (0.5.6)
- Korjattu SAM-tyokalussa ollut virhe (samout generoitiin vasta sitä testissä käyttävän if-lauseen jälkeen).


Muutokset edellisestä versiosta (0.5.5)
- Versiomuutos johtuu Linux-Windows konversion ongelmista. Skripteissä ei muutoksia.


Muutokset edellisestä versiosta (0.5.4)
- norm-affy -skriptissä vsn-normalisointi antoi ulos expressionSetin eikä matriisia. Korjattu skripti vastaamaan
  tätä muutosta R:ssä.
- stat-geneset -skripti ei toiminut current-asetuksella, koska yksi output-tiedosto puuttui. Korjattu.
- stat-linear-model kaatui, jos p-arvojen adjustoinnin käänsi pois päältä. Tämä johtui siitä, että toptable:sta
  yritettiin ottaa sarake P.Values eikä P.Value, joka olisi ollut oikea. Korjattu.
- norm-affy -skriptissä custom.chiptype ei toiminut esimerkiksi RMA-normalisoinnin kanssa. Syynä oli, että
  funktion rma kutsu yritettiin ajaa objektille dat2, jota ei ollut olemassa, kun se olisi pitänyt ajaa
  raakadatan sisältävälle dat-objektille. Korjattu.
- import-soft2 -skriptiin tehty pari päivitystä. Nyt voi valita haluaako datan muuntaa log2-muotoon vain ei.
  Lisäksi korjattu virhe datan kirjoittamisen kanssa. Nyt ekspressioarvot otetaan eset-objektista ulos
  komennolla exprs(eset) eikä suoraan tuota slottia osoittamalla.
- uusi skripti: cluster-kmean-testK, joka testaa klusterointia eri klusterimäärillä, ja antaa tulokseksi
  kuvan, jonka perusteella optimaalista määrää voi arvioida.
- skripteihin norm-agilent, norm-agilent-one-color ja norm-cdna on lisätty normalize.genes -asetukseen myös
  mahdollisuus käyttää VSN-normalisointia.
- NUSE-plottiin lisää tilaa alamarginaaliin, sama heatmap:iin. Lisätty.
- Käytä sample-funktiota HC-kuvan tekemisessä, jos objektien määrä on kovin suuri. Muutettu skriptiä vastaavasti.
- Annotations as a tsv-table. Annotointiskripti palauttaa nyt myös tekstitiedoston, jossa annotaatiot ovat.
- uusi skripti: stat-hyperG-safe, joka tekee GSEA-tyyppisen analyysin KEGG/GO-kategorioille.
- Lisätty Agilent drosophila -annotaatiot normalisointiskriptien valintoihin.
- Jos filtteröintiskripti ajettiin tilastollisen testin jälkeen, eivät p-arvot säilyneet tulostiedostossa.
  Nyt tilanteen pitäisi olla parempi sd/IQR/cv-filtterien osalta, joita ongelma koski.
- Tilastollisiin testeihin lisätty enemmän p-arvojen adjustointiin sopivia menetelmiä (stat-linear-modeling, 
  stat-one-group, stat-several-groups, stat-two-groups). Korjattu myös tähän muutokseen liittyneet virheet
  (tolower muutti kaikki pienille kirjaimille, mikä ei enää toiminut - korjattu if-lauseessa).
- plot-heatmap -skriptissä punainen tarkoitti aliekspressoitunutta ja vihreä yliekspressoitunutta. Vaihdettu 
  toisin päin, mikä on normaalikäytännön mukainen esitystapa.


Muutokset edellisestä versiosta (0.5.3)
- stat-hypergeo -skriptiin kovakoodattu, että KEGG- ja PFAM-paketit ladataan ennen analyysiä. Ennen tämä
  tehtiin automaattisesti, jos analyysiä oltiin tekemässä, mutta Bioconductor-funktio lienee siltä osin
  muuttunut.
- Promoottorianalyysiskriptien polut muutettu Murskan polkuja vastaaviksi.
- Sekvenssianalyysiskriptein polut muutettu Murskan polkuja vastaaviksi.
- stat-snp-chisq -skriptistä puuttui testattavan sarakkeen valintamahdollisuus. Lisätty.
- uusi skripti: promoter-tfbs-cosmo, joka toteuttaa TFBS-haun Bioconductor-projektin Cosmo-pakettia käyttäen.
- Skriptiin seqanal-msa on lisätty mahdollisuus käyttää Clustalia analyysin tekoon, jos dna-paketti on
  asennettu R:ään.
- uusi skripti: norm-illumina-lumi, joka tekee kaikki vaiheet lumi-pakettia käyttäen.
- Korjattu polkumääreet oikein (yksi kauttaviiva puuttui) promoottoriskripteihin.
- Korjattu seqanal-mafft -skriptiin mafft toimimaan (sh-tulkki UNIX:ssa ei tunnista module load komentoa).
- Korjattu seqanal-phylogenetics -skriptiin polkumääre clustaliin oikein (paste-komennosta puuttui pilkku).
- Korjailtu noita seqanal-skriptejä lisää. Nyt pitäisi toimia myös varapalvelimella:lla.
- Korjattu change-interpretation-skriptiin oletusarvojoukkoon log2-linear.
- uusi skripti: delete-columns, joka deletoi datasta valitut sarakkeet. 


Muutokset edellisestä versiosta (0.5.2)
- norm-affy-exon-iptiin lisätty mahdollisuus tuottaa joko geenikohtaiset tai eksonikohtaiset ekspressioarvot.
- norm-affy -skriptiin päivitetty uudet alt CDF:t, joille löytyy suoraan myös annotaatiot osoitteesta:
  http://nugo-r.bioinformatics.nl/NuGO_R.html
- stat-two-group -skriptiin lisätty mahdollisuus verrata kahden populaation variansseja (F-testi)
- stat-linear-model -skriptiin lisätty mahdollisuus olla käyttämättä korjattuja p-arvoja
- simpleaffy-kirjaston kontrolliprobien kuvaukseen on lisätty drosophila2cdf:ää vastaava kuvaus
  (tämä pitänee muuttaa joka kerta kun R päivitetään). Kuvaus löytyy tiedostosta $R_HOME/library/simpleaffy/extdata.
- Weeder-analyysin käyttöön lisätty Drosophila
- promoter-tfbs -skriptiin korjattu virhe (promoottorin pituus: parametreissä short, skriptissä small)
- norm-agilent -skripti kaatuu, jos ajaa sen vain yhdellesirulle. Korjattu (A:sta ja M:stä lopuksi data frame).
   - Sama korjaus tehty skripteihin norm-agilent-one-color ja norm-cdna.
- norm-agilent (ja norm-cdna) -skripteissä spottireplikaattien keskiarvoistaminen aiheutti sen, että rivinimet
  eivät enää tulostuneet, jos normalisoitiin vain yksi siru. Ongelma on nyt korjattu.
- Promoottorityökaluihin lisätty mahdollisuus tehdä analyysi myös Drosophilalle.



Muutokset edellisestä versiosta (0.5.1)
- export-soft -skriptiin korjattu virhe, joka kaatoi export-työkalun, jos datassa ei ollut flagejä.
- norm-affy -skriptiin korjattu pari virhettä. Ensinnäkin MAS%- ja Plier-normalisoinnin jälkeen ei
  palautettukaan log2-muunnettuja arvoja. Toiseksi Plier normalisointi pysähtyi aina virheilmoitukseen,
  sillä custom.chiptype-testissä oli == sen sijaan että siinä olisi oikein lukunut !=.
- norm-agilent-one-color ei toiminut vaan kaatui. Korjattu useita virheitä, jotka johtuivat annotaatioiksi
  merkityn ControlType-sarakkeen väärästä käsittelystä.
   - Samat korjaukset tehty myös norm-agilent-skriptiin.
- norm-agilent-skriptin kontrolliprobien poisto aiheutti sen, että gene set test ei toiminut. Nyt kontrolli-
  probien poisto on vapaaehtoista (default: no), minkä pitäisi korjata tilanne. Replikaattiprobeista otetaan
  silti edelleen keskiarvo.
- uusi skripti: norm-affy-snp, joka ottaa SNP-sirujen CEL-fileet ja laskee niistä genotyypit
- uusi skripti: stat-chisq-snp, joka ottaa normalisoidu SNP-chip datan ja ajaa sille assosiaatioanalyysit
  käyttäen Khiin neliötestiä. Ei ota huomioon perherakenteita.
- filter-by-column -skriptissä on nyt mahdollista testata myös yhtäsuuruus (ennen vain suurempi / pienempi kuin).
- uusi skripti: change-interpretation, joka joko muuntaa log2-arvot alkuperäisiksi tai muuntaa alkuperäiset arvot
  log2-muunnetuiksi.


Muutokset edellisestä versiosta (0.5.0)
- Kaikki tässä versiossa ovat skriptit on tarkoitettu R:n versioon 2.6.1, ja muutokset koskevat vain näitä
  skriptejä!
- Stilisoitu skriptien otsikkohelppejä ja koodien kommentteja.
- Modifioitu seqanal- ja promoter-skriptien polkumäärityksiä siten, että oikea polku tarvitsee osoittaa
  skriptissä vain yhdessä paikassa, joten muuttaminen on helppoa.
- norm-affy -skriptiä korjattu kovalla kädellä. Useimmat normalisoinnit uudelleen annotoiduilla siruilla eivät
  toimineet. Nyt pitäisi toimia.
- plot-venn -skriptissä valittujen tiedostojen nimet eivät tulleet näkyviin. Asia on nyt korjattu ja
  diagrammin eri osat on nimetty tiedostonnimien perusteella.
- calculate-descriptive-statistics -skriptiin lisätty mahdollisuus laskea ne joko geeneille tai siruille.
- norm-illumina -skriptissä oli ongelma: flagejä ei saatu konvertoitua kirjaimiksi. Korjattu ja nyt toimii.
- sort-samples skriptissä oli virhe, eikä se toiminut, jos datassa ei ollut flagejä. Korjattu.
- norm-illumina -skriptiin korjattu Mouse V1-1 paketin nimi oikein.


Muutokset edellisestä versiosta (0.4.9)
- uusi skripti: seqanal-msa, joka tekee usean sekvenssin rinnastuksen FastA-muotoisille sekvensseille.
- uusi skripti: seqanal-phylogenetics, joka tekee fylogeneettisen analyysin Clustal-muotoiselle rinnastukselle.
- Luotu uusi kategoria quality control, johon QC-skriptit on siirretty.
- Fiksattu bugi norm-lme-skriptiin: sanity check oli ennen datan eristystä, ei hyvä. 
- norm-nmds -skriptissä oli kovakoodattu kuvan koko. Muutettu.
- plot-volcano -skriptissä SE-kuvan y-akselin otsikko muutettu muotoon -log(p) eli oikein (oli SE).
- stat-sam -skripti kaatui, jos valitulla Delta-tasolla ei löytynyt yhtään merkitsevää geeniä.
- stat-timeseries -skriptiin tehty muutamia bugikorjauksia (yksi } poistettu, p->p.cut, jne.).
- skriptit, joiden toimivuus R 2.6.1:llä on testattu, ovat tools-R261 kansiossa. Plussalla 
  merkittyihin skripteihin tyhty muutoksia.
   - annotate-genelist
   - average-replicates
   - calculate-descriptive-statistics
   - calculate-fold-change
   - classification-knn
   + cluster-hierarchical
   - cluster-kmeans
   + cluster-qt
   + cluster-som
   - estimate-samplesize
   - export-soft
   - export-tab2mage
   - extract-genes-from-clustering
   - extract-genes-from-stattest
   - extract-samples-from-dataset
   - filter-by-column
   - filter-cv
   - filter-expression
   - filter-flags
   - filter-iqr
   - filter-sd
   - generate-phenodata
   - import-soft2
   - impute
   - merge-tables
   - na-omitted
   + norm-affy
   + norm-affy-exon
   - norm-agilent
   - norm-agilent-one-color
   - norm-cdna
   - norm-illumina
   - norm-lme
   - ordination-nmds
   - ordination-pca
   - (pathway-bayesian)
   - (pathway-boolean-bestim)
   - plot-boxplot
   - plot-chrom-pos
   - plot-corrgram
   - plot-dendrogram
   - plot-heatmap
   - plot-histogram
   - plot-idiogram
   - plot-venn-diagram
   - plot-volcano
   - plot-volcano-data-exists
   - promoter-accessory
   - promoter-cbust
   - promoter-retrprom
   - promoter-tfbs
   + qc-affy
   - qc-affy-rle-nuse
   - qc-agilent
   - qc-agilent-one-color
   - qc-cdna
   - qc-illumina
   - search-correlation
   - search-queryword
   - seqanal-msa
   - seqanal-phylogenetics
   - sort-samples
   - stat-correlate-with-phenodata
   - stat-geneset
   - stat-hyperG-GO
   - stat-hyperG-KEGG-PFAM
   + stat-linear-modelling
   - stat-one-group
   - stat-sam
   - stat-several-groups
   - stat-singleslide
   - stat-timeseries
   + stat-two-groups


Muutokset edellisestä versiosta (0.4.8)
- annotate-genelist -skriptiin oli jäänyt muuttamatta ekspressio- ja p-arvo-sarakkeiden default-arvot
  isoilla kirjaimilla. Tämä korjattu, ja skripti toimii taas.
- affy-norm-exon -skriptillä ei voinut käyttää kuin rma-normalisointia, joten kaikki muut optiot on
  poistettu.


Muutokset edellisestä versiosta (0.4.7)
- uusi skripti: calculate-descriptive-statistics, joka laskee kullekin geenille perustilastoarvot
- uusi skripti: filter-by-column, joka filtteröi datan yhden valitun datan sarakkeen perusteella
- uusi skripti: norm-affy-exon, joka normalisoi Affyn exon-siruja
- norm-agilent-, norm-agilent-two-color- ja norm-cdna -skripteihin on lisätty mahdollisuus käyttää
  taustan korjaamiseen normexp-menetelmää, joka on todettu (Bioinformatics 2007) parhaaksi menetelmäksi,
  kun siihen liitetään pieni offset (50), mikä on myöskin lisätty optioksi skripteihin.
- norm-agilent- ja norm-agilent-two-color -skripteihin lisätty zebrafish-sirun käyttömahdollisuus.
- plot-chromloc -skriptin tuottaman kuvan taustaväriksi vaihdettu valkoinen ja puuttuvien havaintojen 
  väriksi vaalean harmaa, joka tosin rasteroituu ikävästi PNG-kuvassa. Tulos on kuitenkin oletusarvoja
  parempi.
- uusi skripti: sort-samples, joka järjestää näytteet ja fenodatan jonkin fenodatan sarakkeen mukaan
- norm-agilent-skriptiä muutettu siten, että se ei enää kaadu, vaikkei ControlType-saraketta olisikaan
  tuotu sisään.
    - sama muutos tehty myös norm-agilent-one-color -skriptiin 


Muutokset edellisestä versiosta (0.4.6)
- geneset test: tulosta tulokset myös listana, jotta näkee kaikkien kategorioiden tulokset. Korjattu.
- norm-agilent (1- ja 2-väri):
   - Replikaattispoteista lasketaan nyt oletusarvoisesti (eikä sitä voi kääntää pois päältä) keskiarvo,
     jos spotit ovat samalla sirulla.
   - Poistetaan kontrollispotit datasta automaattisesti
- norm-cdna:
   - Replikaattispoteista lasketaan nyt oletusarvoisesti (eikä sitä voi kääntää pois päältä) keskiarvo,
     jos spotit ovat samalla sirulla.
- PCA - add a possibility to scale the value before analysis. Lisätty.
- uusi skripti: ordination-nmds, joka tekee siruille NMDS-analyysin, ja piirtää siitä kuvan
- uusi skripti: filter-cv, joka filtteröi geenit niiden CV:n perusteella
- QT- ja K-means -clustering: R:n tuottamissa kuvissa eivät y-akselit ole samassa skaalassa. Korjattu.
- uusi skripti: qc-agilent, joka tekee QC:n Agilentin 2-väridatalle
- uusi skripti: qc-agilent-one-color, joka tekee QC:n Agilentin 1-väridatalle
- korjattu bugi stat-group-testistä. Empirical Bayes -testistä puuttui mahdollisuus bonferroni-korjauksen
  käyttöön (sitä ei laskettu skriptissä ollenkaan, vaikka se oli valittu optioksi).
- norm-agilent -skriptissä oli bugi sirutyyppi mgug4112a:n kohdalla: a puuttui lopusta. Korjattu.
- calculate-fold-change -skriptissä ollut fold change laskenta korjattu oikein (FD<-dat3[,2]/dat3[,1]
  muutettu muotoon FD<-dat3[,2]-dat3[,1]).


Muutokset edellisestä versiosta (0.4.5)
- empirical bayes -testien pitäisi tulostaa myös fold change (ei pelkästään lineaarisen mallin)
  Korjattu kahden ryhmän testin osalta. Usean ryhmän tapauksessa pitäisi käyttää linear modelling-
  työkalua.
- limma-analyysi ei tulosta rivinimiä -> korjattava ASAP! Korjattu. Lisäksi muutettu design-matriisin
  kirjoitusta siten, että levylle kirjoitettava tiedosto tulkitaan pelkkänä tekstinä.
- limma-analyysi laittoi geenit väärään järjestykseen (tulokset on sortattu suurimmasta pieneimpään,
  mitä ei ollut otettu huomioon tuloksia levylle kirjoitettaessa). Korjattu, toivottavasti oikein.
- annotate-skriptiä päivitetty siten, että jos datalle on tehty tilastollisia testejä siten, että p-arvo
  ja fold change ovat saatavissa, niin käyttäjä voi valita lisäävänsä ne annotoituun geenilistaan
  lisäinformaatioksi.


Muutokset edellisestä versiosta (0.4.4)
- Illumina normalisointi ei toiminut, koska sarakkeiden nimet datassa ja rivien nimet fenodatassa
  eivät sisältäneet tiedostopäätettä. Nyt asian pitäisi olla kunnossa, ja homman toimia.
- stat-two-groups kaatui, koska sanity check-vaihe oli ennen muuttujien luomista, joiden suhteen
  nuo checkit ajettiin. Korjattu.
- Extract-samples-from-dataset kaatui, jso datasetti ei sisältänyt call-arvoja. Nyt pitäisi toimia
  myös datalle, jossa call-arvoja ei ole.
- stat-several-groups kaatui Kruskall-Wallisin testiin, sillä muuttujaa p ei ollut alustettu ennen
  silmukkaa. Korjattu.
- QT-clustering kuva R:stä saatuna näyttää pahalta (x-akseli on liian pitkä, fiksaa kuten K-means:ssa).
  Korjattu.
- correlate-with-phenodata kaatuu, jos datasta on ensin laskettu ryhmäkeskiarvot
   - tsekkaa, että datan ja fenodatan dimensiot ovat samat - jolleivat, anna virheilmoitus.
   - Korjattu.
- cluster-som kaatui siihen, että dists-objektista otettiin ensimmäinen sarake, vaikka objekti oli 
  oikeasti vektori. Korjattu.


Muutokset edellisestä versiosta (0.4.3)
- calculate-fold-change -skriptiä muutettu hieman. Nyt se ei laskenut fold changea edes kahdelle ryhmälle.
  Lisäksi ulostulotiedoston nimeksi muutettu VVSADL:n vaatima fold-change. Lisäksi fold change talletetaan
  sarakkeeseen chip.FD entisen FD:n sijaan.
- plot-heatmap -skriptiä korjattu siten, että se odottaa tulostiedostoksi heatmap.png eikä esim.
  dendrogram-color.png, kuten aiemmin.
- plot-dendrogram -skriptiin lisätty phenodata, jota ilman skripti ei tuottanut dendrogrammia.
- Muutettu tilastoskriptien oletusarvoiseksi metadata-sarakkeeksi group.
- stat-correlate-with-phenodata -skriptistä korjattu virhe: yritettiin pyöristää myös flagejä
  korrelaatioarvojen lisäksi. Skripti kaatui virheilmoitukseen, että merkkidataa ei voi pyöristää.


Muutokset edellisestä versiosta (0.4.2)
- Kaikki merge-funktiot poistettu yhtä geneeristä lukuunottamatta.
- norm-illumina -skripti muutettu toimimaan siten, että data on tuotu
  syöttötyökalulla, ja se on kirjoitettu levylle yksittäisiksi tiedostoiksi.
- stat-two-groups- ja stat-several-groups -skripteistä poistettu viittaus SAM-
  analyysiin. 
- uusi skripti: stat-sam, joka laskee SAM-analyysi yhdelle, kahdelle tai useammalle
  ryhmälle.
- KNN ei toimi, sillä group-saraketta ei voi jättää tyhjäksi test datalle.
   - virhe johtuu siitä, että predict-funktiolle syötetään train-setti, joka
     on oikean mittainen, mutta cl-muuttuja puolestaan on koko vektorin mittainen
     (=ei hjuva).
   - Korjattu.
- filter-IQR ei tuota tulostiedostoa (library-komento puuttuu). Korjattu.
- K-means kuva, joka tuotettu R:ssä näyttää oudolta. Nyt näyttää paremmalta. Koodi
  oli skriptissä hieman väärässä järjestyksessä, mutta nyt pitäisi olla oikein (
  kuvan piirto ennen taulukon kirjoitusta).
- calculate-fold-change kaatuu. Ongelmana ainakin vaadittava ryhmien määrä. Nyt 
  pitäisi toimia kunnolla.
- plot-dendrogram kaatui, koska toisen ulostulotiedoston nimi oli väärin.
- plot-corrgram kaatui, koska corrgram-kirjastoa ei ollut ladattu ennen analyysiä
- plot-heatmap kaatui, koska heatmap() vaati syötteksi matriisin, ja sille annettiin
  data frame. Korjattu.
- Jos niissä skripteissä, joissa pitää valita jokin fenodatan sarake ei ole valittu
  mitään saraketta, niin ajetaan skripti group-sarakkeella, joka pitäisi aina olla
  täytettynä. Muutokset tehty skripteihin:
     - average-replicates
     - calculate-fold-change
     - extract-samples-from-dataset
     - extract-genes-from-stattest
     - norm-lme
     - plot-dendrogram
     - plot-volcano-data-exists
     - stat-correlate-with-phenodata
     - stat-geneset
     - stat-linear-modelling
     - stat-sam
     - stat-several-groups
     - stat-timeseries
     - stat-two-groups
     

Muutokset edellisestä versiosta (0.4.1)
- stat-linear-modelling -skripti kirjoittaa nyt ulos myös kunkin parametrin suhteen
  lasketun fold changen.
- uusi skripti: extract-genes-from-stattest, joka erottaa tilastollisten testien
  antamista tuloksista ne, joiden p-arvo on alle valitun raja-arvon.
- uusi skripti: volcano-plot-data-exists, joka plottaa tehdyn tilastollisen testin
  (linear model tai empirical bayes) tulokset Volcano plotin muotoon.
- Muutettu muutamiin plottiskripteihin kuvan koon asetukset oikein (600/72 -> w/72 ym.).
- norm-illumina -skriptiin lisätty mahdollisuus käyttää myös ProbeID-annotaatioita.
- norm-illumina -skriptistä korjattu flagging-muuttujan käsittelyyn liittyviä virheitä.


Muutokset edellisestä versiosta (0.4.0)
- uusi skripti: stat-linear-modelling, jossa on toteutettu lineaaristen mallien käyttö data-
  analyysiin, kuten se on limma-paketissa toteutettu.
- norm-illumina-skriptiin lisätty mahdollisuus valita sarake-erotin (, tai ; tai tab). Jos
  erotin on ;, on pilkun merkkinä , eli käytetään "suomalaisia asetuksia".
- uusi skripti: norm-agilent-average-replicates, joka laskee normalisoinnin jälkeen keskiarvon
  sirulla olevista kahdesta täpläreplikaatista. Tämä toimii vain, jos replikaatteja on tismalleen
  kaksi.
- stat-skripteihin lisätty parametri, jolla käyttäjä voi valita, minkä fenodatan sarakkeen
  suhteen testi halutaan tehdä.
- Hypergeometrinen testi ei palauttanut ikinä yhtään merkitsevää tulosta, koska skriptille
  annetun datan sijaan testeissa käytettiin sirun kaikkien geenien listaa -> kaikilla geeneillä
  p-arvo 1. Korjattu.


Muutokset edellisestä versiosta (0.3.9)
- norm-lme-skripti muutettu siten, että käyttäjä voi valita random-efekti sarakkeen.
- cluster-som -skripti tuottaa nyt tuloksista myös R:ssä generoidun kuvan.
- cluster-kmeans -skripti tuottaa nyt myös kuvan klustereista.
- cluster-qt -skripti tuottaa nyt myös kuvan klustereista.
- uusi skripti: extract-genes-from-clustering, joka tallentaa uudeksi datatiedostoksi
  tietyn klusterin geenit.
- uusi skripti: norm-agilent-one-color, joka normalisoi Agilentin 1-väridatan.
  (elina.palonen@btk.fi).
- uusi skripti: norm-illumina-separate-files, joka normalisoi Illumina sirut, jos ne
  on tuotu sisään import tool:n kautta.
- uusi skripti: plot-corrgram, joka visualisoi näytteiden väliset korrelaatiot hauskana
  kuvana. 
- uusi-skripti: extract-samples-from-dataset, joka antaa käyttäjän poistaa näytteitä ja
  fenodatan rivejä tietyn fenodatan sarakkeen perusteella.
- normalisoimaton Illumina-data tuli ulos raaka-arvoina. Muutettu tuottamaan log2-muunnettuja
  arvoja.
- uusi skripti: merge-tables-exp-with-annot, joka yhdistää tabulaarisen ekspressiodatan
  tabulaarisen annotaatiodatan kanssa.
- norm-illumina -skriptiä päivitetty siten, että ainoastaan jos flagit halutaan tuottaa,
  oletetaan että datassa on olemassa Detection-sarakkeet. Lisäksi, jos Detection sarakkeita ei
  ole, ei skripti edes yritä tuottaa flagejä. Pitäisi olla vikasietoisempi nyt.
- uusi skripti: filter-iqr, joka filtteröi geenit niiden kvartiilivälin pituuden perusteella.


Muutokset edellisestä versiosta (0.3.8)
- stat-hyperG -skripti jaettu kahtia siten, että GO-analyysi on omassa skriptissään, ja
  KEGG/PFAM omassaan.
- norm-agilent -skripti ei kirjoittanut flag-arvoja, vaikka ne oli dataa tuotaessa määritelty.
  Tämä muutettu toimimaan samalla periaatteella kuin cDNA-normalisoinnissa.
- filter-sd ei toimi yhdellä sirulla, eikä toimi filter-flag:kään. Tällöin yritetään
  valita subscriptillä tyyliin [1,] vektorista tai tyyliin [,1] vektorista. Paha, paha juttu.
   - Nyt tämä on korjattu, ja skriptien pitäisi toimia. Lisäksi nyt voi SD-filtterissä käyttää
     vain geenien suhteen filtteröintiä, joten siruja on oltava vähintään 2.
- uusi skripti: na-omit.R. Poistaa datasta kaikki rivit, joilla on puuttuvia havaintoja.
- cDNA- ja Agilent-normalisointi ei nyt enää kaadu flag-sarakkeen puuttumiseen, mutta 
  kylläkin, jos jokin muista vaadituista sarakkeista puuttuu.
- Käyttöliittymä kaatuu, jos stat-hyperG-GO ajetaan, eikä yhtään merkitsevää geeniä löytynyt. 
  Skriptiä muutettu hieman toisenlaiseksi, jolloin HTML-file luodaan suoraan kirjoittamalla
  tagit tiedostoon eikä R:n helpperi-funktiolla.
- SD-filtterissä oli edelleen maininta muuttujasta marg. Se on nyt poistettu if-lauseesta,
  ja nyt skriptin pitäisi toimia.
- uusi skripti: plot-venn-diagram piirtää Vennin diagrammin kolmelle datasetille.
- uusi skripti: calculate-fold-change.R laskee kahden ryhmän keskiarvojen perusteella fold changen
- PCA-skriptissä ollut virhe korjattu (samples->chips). Nyt PCA toimii myös näytteille.
- uusi skripti: merge-illumina-files, joka yhdistää kaksi systeemiin ladattua bead summary data 
  -tiedostoa. Yhdistämisen jälkeen tiedoston voi normalisoida norm-illumina -skriptillä.
- uusi skripti: generate-phenodata, joka kirjoittaa tyhjän phenodatan, jos skripti ajetaan
  esinormalisoidulle datalle.
- uusi skripti: plot-boxplot, joka piirtää boxplotin normalisoidulle datalle.
- stat-geneset -skriptissä viitattiin vielä hgu133aSYMBOL-ympäristöön, joka korvattiin siten,
  että skripti toimii kaikilla annotaatiopaketeilla.
- search-queryword-skriptissä viitattiin edelleen hgu133a-ympäristöön, mutta tämä muutettiin käyttämään
  oikeaa ympäristöä kaikille siruille, joista annotaatioita on.
- Kaikkin plot-funktioihin lisätty mahdollisuus säätää kuvan kokoa.
- norm-cdna- ja norm-agilent -skripteihin tehty pieni muutos flagien käsittelyyn, jotta ne saadaan mukaan
  aina, olipa sarakkeet määritelty miten vain. 
- uusi skripti: plot-dendrogram, joka piirtää kaksi erilaista dendrogrammia.
- uusi skripti: plot-heatmap, joka piirtää heatmapin.
- uusi skripti: plot-histogram, joka piirtää kullekin sirulle histogrammin


Muutokset edellisestä versiosta (0.3.7)
- Manuaalia laajennettu ja parannettu huomattavasti.
- Illumina RatRef-12 lisätty annotaatiokirjastoihin, ja normalisointiskriptiä
  korjattu vastaavasti.
- Agilent-sirujen (ihminen, hiiri ja rotta) annotaatiot lisätty, ja normalisointiskriptiä
  muutettu vastaavasti. Nyt Agilent-normalisoinnin yhteydessä pitää valita sirutyyppi.
- cDNA-QC-skriptin kaikkien MA-plottien Y-akselien pitäisi olla samalla skaalalla. OK.
- Käännetty qc-affy-rle-nuse -skriptin sirujen nimet pystysuoraan, jotta ne ovat luettavissa
  kaikille siruille.
- ClusterBuster-skripti palauttaa nyt tuloksessa transkriptiofaktoreitten nimet, ei niiden
  matriisiluokkia kuten aiemmin.
- Muuta fenodata tuottamaan sarakkeet järjestyksessä: sample, chiptype, group, training. OK.
- hyperG-skriptin silmukoissa oli pientä säätämistä, koska if..else-lauseet oli muotoiltu väärin.
  Nyt ne on korjattu, ja skriptin pitäisi toimia oikein.
- Volcano plot implementoitu.
- Agilent-normalisointia korjattu - sirutyypin valinta ei onnistunut. Nyt pitäisi toimia.


Muutokset edellisestä versiosta (0.3.6)
- Lisätty skriptiin stat-several-groups parametrin use.simple.analysis, jota AffyWizard
  käyttää asetuksella yes. Oletusarvoisesti se on pois päältä.
    - stat-several-groups-foAffyWizard skripti on poistettu.
- importGDS-skripti lisää nyt sarakkeiden nimien eteen "chip.", jotta analyysitkin toimivat.
- stat-geneset -skriptistä korjattu pieni bugi, joka esti kuvan otsikoiden tulostumisen oikein
  oli (names(substr(...)) kun pitäisi olla subtr(names(...)).
- stat-hyperG -skriptistä poistettu yksi ylimääräinen } yhden if-lauseen lopusta.
- search-correlation -skriptistä poistettu debuggaus-rivit.
- Päivitetty hc-skriptiä: write.tree() -> format="Newick" poistettu.
  - Muutettu takaisin, sillä muutoin skripti antaa jotain ihan höpöä.
- Päivitetty skriptiä stat-hyperg siten, että sen pitäisi toimia, vaikkei merkitseviä
  tuloksia olekaan saatu.
- Lisätty norm-cdna -skriptin generoiman fenodatan sample-sarakkeen nimien perään ".tsv". Nyt
  normalisoinnin pitäisi toimia.
- suurenna hc-bootstrapping fonttikokoa - OK, nyt 0.75 (ennen 0.5)
- PCA kaatuu (Error in if (no < 3) { : missing value where TRUE/FALSE needed)
    - Ongelmana oli parametri, joka pyydettiin prosentteina (0-100) eikä desimaalilukuna
      (0-1). Ongelma korjattu jakamalla prosenttiluku sadalla.
- average-replicates -skriptin pitää palauttaa sarakeotsikot tyyliin "chip.group1" ei
  "group1" - OK.
- norm-cdna -skripti ei kirjoita flagejä, vaikka ne olisi määritettykin - OK
   - Flagit luettiin väärästä objektista (dat3), mutta nyt sen pitäisi toimia (luetaan
     objektista dat)
- impute-skriptiin lisätty tarkistus, että onko A (average) -matriisin sarakeluku yli 0.
  Jos on, niin kirjotettavaan tiedostoon liitetään myös average-sarakkeet, joten imputen
  pitäisi toimia nyt myös cDNA-datalle.
- qc-cdna kaatuu (mihin?)
   - Useita juttuja meni pieleen, erityisesti puuttuvien havaintojen tapauksessa.
     Nyt pitäisi toimia.
- filter-expression: korjattu pieni bugi, joka esti funktion toimimisen yhdellä sirulla
  (tällöin palautettiin scaled.dat2-vektori, jolle apply-funktiota ei voinut käyttää, nyt
  muutetaan tuo vektori ensin dataframeksi).
   - sama ongelma vaivasi filter-sd -skriptiä, mutta sekin on nyt korjattu.
- filter-flags antoi virheilmoituksen yhdelle sirulle, koska vektorin sarakkeiden lukumäärä
  oli looginen nolla. Nyt muutettu flagit ensin datakehikoksi, ja sitten siitä lasketaan
  sarakkeiden lukumäärä.


Muutokset edellisestä versiosta (0.3.5)
- Lisätty stat-several-groups- skriptiin empirical Bayes -menetelmä (limma).
  Affy-wizard käyttää tätä, ja se toimii todella nopeasti.
- stat-two-groups -skriptiin lisätty empirical Bayes -menetelmä (limma)
- norm-cdna-skripti luki rivien nimet annotaatioksi. Muutettu rivien nimet geenien nimiksi
  ja annotation-sarake annotaatioksi.
- Utility-kategoriaan väärin sijoitetut skriptit on nyt siirretty utilities-kategoriaan.
- ClusterBuster palautti tyhjän tiedoston, koska jasparin matriisitietokantaa ei löytynyt. 
  Ongelma on nyt korjattu (jaspar_cor.tsv muutettu ohjelmakutsussa oikein muotoon 
  jaspar_core.txt).
- export-funktiot kaatuivat heti alkuun, koska luettavan fileen nimi oli normalized. Oikea
  file oli tietysti normalized.tsv, mikä onkin nyt korjattu.
- Lisätty cbust-skriptiin ohjelman asetukset (-c, -m, ... ,-p)
- Lisätty pudotusvalikko norm-illumina-skriptiin. Valikosta valitaan siru, jonka mukaan valitaan
  käytettävä annotaatiokirjasto.
- search-queryword -skriptistä korjattu kaksi bugia, joiden takia haku aina epäonnistui.
  Query oli muodoissa "chr" ja "genename". Nyt yhtenäistetty.
- Affy-wizardille tehty oma versio stat-several-groups -skriptistä. Se palauttaa merkitsevät
  geenit, jos niiden adjusted p-arvo on yli p-arvon kriittisen arvon. Muutoin palautetaan
  100 raa'alta p-arvoltaan merkitsevintä geeniä.
- Muutettu one-sample-text raportoimaan adj.p.value -sarake entisen sekasikiönimen sijaan.
- Muutettu estimate-samplesize tuottamaan 600/72 tuuman kokoisia kuvia entisen 600 tuuman sijaan.
- clust-hc -skriptissä ollut kirjastojen latausvirhe korjattu, ja nyt bootstrap-analyysinkin
  pitäisi toimia.
- Kommentoitu SOM-skriptistä pois gridin koon järkevyystestaus.
- stat-genesettest-skriptissä poistettu parametri ylab par-asetuksesta ja siirretty se plot-komentoon,
  mihin se oikeasti kuuluisikin.


Muutokset edellisestä versiosta (0.3.4)
- Muutettu .txt -> .tsv, joka viittaa nyt vain tab-delimited text -tiedostoon. Samaten
  muutettu .text -> .txt, joka viittaa nyt geneeriseen tekstitiedostoon.
- Muutettu SOM-skriptin oletusväritys sini-keltaiseksi
- Muutettu hierarkkisen ryhmittelyanalyysin skriptiä siten, että tuloksena
  on bootstrapin tapauksessa 600/72 tuuman kuva 600 tuuman kuvan sijaan.
- Lisätty parametrien tsekkausta skripteihin:
   - average-replicates
   - classification-knn
   - cluster-hierarchical
   - cluster-kmeans
   - cluster-som
   - filter-flags
   - filter-sd
   - norm-lme
   - stat-one-group
   - stat-several-groups
   - stat-two-groups   
- export-soft -funktio (kuten tab2mage-funktiokin) tulostaa nyt vain blankon SOFT-fileen.
- fiksattu bugi impute -skriptissä
- fiksattu bugi stat-correlate-with.phenodata -skriptissä. Ajo kaatui, koska skriptille ei
  syötetty fenodataa, vaikka se tarvittiin skriptin suoritukseen.
- Hierarkkinen clusterointi palauttaa nyt tiedoston hc.txt entisen hc.tsv:n sijaan.
- norm-affy -skriptissä valitaan nyt custom chip type pudotusvalikosta vahinkojen välttämiseksi.


Muutokset edellisestä versiosta (0.3.3)
- Normalisointiskripteihin tehty muutos siten, että skripti kirjoittaa fenodatan 
  käyttöliittymän sijaan. Samalla lisätty parametri, jolla voidaan asettaa custom
  chiptype. Muutokset tehty skripteihin:
    - norm-affy.R
    - norm-illumina.R
    - norm-cdna.R
    - norm-agilent.R
- Normalisointiskriptit tuottavat nyt tiedoston phenodata.txt eivätkö yritä enää
  lukea sitä sisään, ja tuottaa phenodata_out.txt:tä.
- Uusi skripti: tab2mage.R, joka kirjoitta tyhjän tab2mage-tiedoston, jonka käyttäjä
  voi täyttää vaikka Excelissä.
- Uusi skripti: import-SOFT.R, joka hakee GDS:n suoraan GEO:sta ko. numerolla, ja
  luo siitä normalisoidun datan ja fenodatan.
- Illumina koodia muokattu siten, että se pystyy lukemaan myös BeadStudio 3:lla tuotetut 
  tiedostot, joissa header-informaatiota on 8 rivia aiempien versioiden 7 sijaan. Tämä
  vaati pieniä muutoksia datan sisäänluvussa ja flag-arvojen generoinnissa. Lisäksi
  käyttäjän on nyt itse valittava BeadStudion versionumero, jolla data on generoitu.


Muutokset edellisestä versiosta (0.3.2)
- Uusi skripti: promoter-retrprom, joka vain hakee promoottorisekvenssit tietokannasta,
  muttei tee niillä se pidemmälle meneviä analyysejä.
- qc-affy -skripti muokattu käyttämään simpleaffy-kirjaston funktioita.
- Uusi skripti: qc-affy-rle-nuse, joka tuottaa RLE- ja NUSE- plotit Affymetrix datalle.
- Uusi skripti: stat-correlate-with-phenodata, joka korreloi geeniekspression 
  phenodatan kanssa. 
- Uusi skripti: average-replicates, joka laskee kunkin ryhmän keskiarvon.
- Uusi skripti: stat-hyperG, joka tekee hypergeometrisen testin GO/KEGG/PFAM-luokille.


Muutokset edellisestä versiosta (0.3.1)
- qc-affy -skripti piirtää nyt eri sirut eri väreillä, ja tekee kuvaan kuvatekstin,
  jossa sirujen nimet ja värit on spesifioitu.
- clust-som -skrptiin lisätty sarake, joka kertoo gridin koon (x*y). GUI:n piirtämä
  kuva tehdään tämän perusteella.
- Korjattu promoter-cbust -skriptissä ollut virhe, joka esti promoottorisekvenssejä
  latautumasta oikein objektiin upstream.
- Korjattu search-correlation -skriptissä ollut vika, joka kaatoi analyysin, jos 
  grep-komento löysi useampia geenejä, joissa sama string esiintyi.
- Korjattu search-queryword -skriptissä parametrien nimiin liittyvä ristiriita.
  Näin ainakin AFfymetrixID:n perusteella toimiva etsintä saatiin toimimaan.
- promoter-cbust -skriptissä ollut tulostiedosto clusters.txt muutettu muotoon
  clusters.text, jolloin käyttöliittymä osaa sen näyttää (.txt viittaa aina matriisiin).


Muutokset edellisestä versiosta (0.3.0)
- norm-illum -skripti muutettu toimimaan käytetyn Bioconductor-version (1.9) kanssa.
- norm-illum -skriptin normalisointifunktiot (qspline ja quantile) muutettu käyttämään 
  affy-paketin normalisointia.
- lisätty skripti Agilent-datan sisäänlukua ja normalisointia varten.
- cDNA / Agilent -skripteihin tehty muutoksia flagien käsittelyn osalta. Jos flagejä ei ole
  ei niitä yritetä enää kirjoittaakaan, joten skripti ei kaadu niihin.
- Korjattu cbust-skriptissä kirjoitusvirhe, joka esti oikean annotaatiokirjaston latautumisen.
- filter-sd -skriptin kuvausta muokattu: retain->filter out.
- stat-geneset -skriptissä muokattu kuvien marginaaleja.


Muutokset edellisestä versiosta (0.2.9)
- Hierarkkinen klusterointi tuotti dendrogram.txt:n, mutta GUI odotti hc.txt:tä. Korjattu.
- KNN-skriptissä parametrien arvojen testauksessa oli kirjoitusvirhe (sulku puuttui). Korjattu.
- Idiogram menee ajoon, ja tekee tuloksen, mutta kuva ei ilmesty. GUI odotti tiedostoa
  idiogram.png, mutta tuotettiinkin midiogram.png, mistä ongelma johtui. Korjattu.
- Geneset test -skriptiin korjattu GO-luokkien testaus toimivaksi (muuttujan nimi oli väärin.
- Korjattu bugi plier-kirjaston latauksessa plier-normalisoinnin yhteydessä.
- KEGG/GO-analyysin tuottaman kuvan marginaaleja pienennetty, jotta kuvasta tulee informatiivisempi.
  Samalla myös tekstikokoa suurennettu hieman.
- qc-affy -skripti korjattu toimimaan myös yhdellä sirulla.
- qc-affy -skriptiä muutettu siten, että sirut ovat nyt sarakkeina ja kontrolliarvot riveinä. 
  Lisäksi muutettu rivien ja sarakkeiden nimiä hieman.
- qc-affy-skriptissä oletettiin tulevan .jpg-kuva, mutta tulikin .png-kuva. Korjattu.
- stat-geneset-skriptistä korjattu kirjoitusvirhe, joka esti kuvien marginaalien skaalauksen.
- SD-filtteriskriptissä "retain" muutettu muotoon "filter.out".
- promoter-cbust -skriptissä korjattu promoottorien kokoon liittyvän vertailuoperaattorin
  kirjoitusvirhe.
- promoter-tfbs -skriptiä muutettu siten, että medium-analyysi ajaa myös small-analyysin. 
- promoter-tfbs -skriptistä poistettu mahdollisuus ajaa large- ja extra-analyysejä
- promoter-tfbs -skripti ei enää kirjoita seqs.txt.wee tiedostoa. Tai kirjoittaa kyllä,
  mutta sitä ei enää näytetä käyttäjälle.


Muutokset edellisestä versiosta (0.2.8)
- KNN-skripti muutettu toimivaksi (erottelee .chip- ja muut sarakkeet toisistaan aluksi)
- Hierarkkinen klusterointi muutettu käyttämään pvclust-kirjastoa bootstrap-analyysissä.
- Hierarkkinen klusterointi tuottaa sulkukaavion, joka voidaan parsia VISKI-projektin
  visualisointipalikan käyttöön.
- Search-correlation-skriptissä korjattu virhe, joka aiheutti sen, että skripti toimi vain
  Affyn hgu133a-siruilla, ja löysi silloinkin vain tietyn geenin kanssa samalla tavalla
  ekspressoituneet geenit.
- Idiogrammi-funktiolta puuttui tieto ihmisen kromosomiraidoista. Ongelma korjattu, ja nyt
  skripti toimii.
- Fiksattu bugi one-group test ja several groups test -skripteissä. Tämä liittyi korjattujen
  p-arvojen laskemiseen.
- Fiksattu pieni bugi (?) PCA-skriptissä. Jostain syystä skripti toimii hyvin R:ssä, muttei enää
  jos se suoritetaan NAMI:n kautta. Testataan toimiiko korjattukaan versio NAMI:ssa.
     - Ei toimi vieläkään, vaikka samalla datalla skripti toimii R:ssä
- Fiksattu SOM-skriptissä värien määrä vastaamaan solujen lukumäärää siten, että kukin solu,
  joka on yhtä kaukana ensimmäisestä solusta värjätään samalla värillä.
- Chromosome position plot:ssa ei eroteltu dataa ja flagejä, joten skripti kaatui heti kärkeen.
- Muutamien kuvien tuotossa tiedostopäätteet muutettu oikeiksi (jpg->png).
- Stat-two-group-tests -skriptissä ladataan nyt LPE-kirjasto ennen ko. testin suorittamista.
- Stat-one-group -skriptissä muutettu Wilcoxon-testin käyttämä datakehikko vektoriksi, jolloin
  testi toimii oikein.
- Search-queryword -skriptissä korjattu verailuoperaattoriin liittyvä virhe, joka esti skriptin 
  ajamisen.
- Muutettu search-correlation -skriptiä siten, että korrelaatio voi olla väliltä 0-1 (aiemmin
  -1...1).
- Kutsut png-komentoon on muutettu muotoon bitmap("filename", width=w/72, height=h/72).
- Korjattu muutama lyöntivirhe (dat, kun piti olla dat2) plot-chrom-pos -skriptissä.
 

Muutokset edellisestä versiosta (0.2.7)
- Manuaalin ensimmäinen draft valmis.
- Muutettu hakemistopuuta hieman
- Muutettu normalisoinnin yhteydessä tapahtuvaa phenodatan levylle kirjoitusta: phenodata-taulukko kirjoitetaan aina
  aluksi levylle, ja siihen lisätään myöhemmin chiptype-sarake, jos se siitä puuttuu.


Muutokset edellisestä versiosta (0.2.6)
- Skriptien nimiä GUI:ssa muokattu hiukan (suodin-skripteissä kerrotaan tämä nyt jo nimessä, jne.)
- Affymetrix-normalisointeihin lisätty PLIER-algoritmi
- MAS5- ja PLIER-normalisoidulle datalle on mahdollista tehdä VSN-normalisointi (mutta muille ei, sillä niiden
  antamat arvot on jo log2-muunnettu, mikä ajaa saman asian)
- Lisätty skriptiin classification-knn.R tarkistuksia asetusten oikeellisuuden suhteen ja poistettu phenodata-taulukkoa
  muokkaavat ominaisuudet.
- Lisätty skripti otoskoon ja voiman arviointiin (estimate-samplesize.R)
- Hierarkkinen klusterointi toteutettu uudelleen amap-kirjaston hcluster-funktiota käyttäen. Tämä mahdollistaa
  puun muodostamisen suurillekin (>20000 havaintoa) aineistoille suhteellisen lyhyessä ajassa (~30 min).
- Hierarkkisen klusteroinnin uudelleenotantafunktio on muutettu täysin. Nykyinen toteutus käyttää pvclust-kirjastoa,
  joka tuottaa valmiin kuvan, johon bootsrapping-arvot on aseteltu. Aiempi toteutus laski bootstrapping-puun hitaasti
  yksittäisiä funktioita käyttäen. Nykytoteutus on nopeampi (yli 6000X), mutta toisaalta ainoastaan pearsonin 
  korrelaatiota voidaan käyttää puun muodostamiseen.
- Tunnettujen transkriptiofaktoreitten sitoutumiskohtien etsimiseen sekvensseistä käytetään ClusterBuster-ohjelmaa.
  Tämä on uusi skripti TFBS-työkalujen joukkoon.


Muutokset edellisestä versiosta (0.2.5)
- GEO:ssa tapahtuneet batch submission -ohjeiden mukaiset muutokset toteutettu export-soft.R skriptiin.
- Idiogrammin piirtävä skripti muokattu toimivaksi
- Geenien kromosomaaliset positiot visualisoiva skripti plottaa nyt eri väreillä yli- ja aliekspressoituneet geenit
- Phenodataan ei enää kirjoiteta chiptype-saraketta, jos se on jo olemassa


Muutokset edellisestä versiosta (0.2.4)
- Promoottorityökalu muokattu toimimaan UCSC:stä ladatuilla RefSeq-promoottorisekvensseillä
   - weederlauncher.out ei toiminut, joten skripti käyttää nyt weederTFBS:ää
- Korjattu bugi filter-expression.R skriptissä. Bugi aiheutti sen, ettei skripti poistanut
  filtteröinnissä datasta yhtään geeniä. 


Muutokset edellisestä versiosta (0.2.3)
- LME-skripti muutettu käyttämään nlme-kirjastoa lme4:n sijaan.
- bestim.R-skripti tehty bnBestFit-työkalun integrointia varten.


Muutokset edellisestä versiosta (0.2.2)
- Uudet skriptit:
   - cPlot-cColor.R
   - GOFisherTest.R
   - idiogram.R
   - promoter.R
- LME-skripti korjattu käyttämään oikeita residuaaleja


Muutokset edellisestä versiosta (0.2.1):
- Skriptien parametrien nimet on pyritty muuttamaan informatiivisiksi


Muutokset edellisestä versiosta (0.2.0):
- Phenodata luetaan samasta kansiosta kuin datakin
- Bayesian network-skripti on omana kategorianaan eikä enää timeseries-skriptien osana
- One-sample t-testi toimii noin 500-kertaa nopeammin kuin aiemmin, funktio kirjoitettu uudelleen


Muutokset edellisestä versiosta (0.1.4):
- Bootstrapping-validointi lisätty hierarkkiseen ryhmitttelyanalyysiin
- Aikasarja-analyysi:
   - Etsi periodisesti muuttuvat geenit
   - Muodosta network
   - Independent component analysis


Muutokset edellisestä versiosta (0.1.3):
- Muutettu otsikkotietoja edelleen...


Muutokset edellisestä versiosta (0.1.2):
- Korjattu otsikkotiedot
   - Kirjoitusvirheet korjattu
   - Skriptit luokiteltu analyysiympäristön luokkiin:
      - Normalisation
      - Preprocessing
      - Clustering
      - Statistics
      - Utility
   - Käyttämättömät luokat
      - Spot filtering
      - Expression
      - Pathways


Muutokset edellisestä versiosta (0.1.1):
- Linear Mixed Model korjattu
   - Mahdollistaa satunnaisvaikutusten (batch effect) poistamisen normalisoidusta datasta


Muutokset edellisestä versiosta (0.0.2):
- tuki cDNA- ja Illumina-siruille
   - normalisointi
   - monet analyysioptiot toimivat myös näillä siruilla
- tuki phenodatalle
   - kaikki skriptit lukevat tarvittaessa parametrit tiedostosta phenodata.txt
   - phenodata.txt:n oletetaan olevan hakemistossa ../../common. Jos se puuttuu
     kirjoitetaan normalisointiskriptien yhteydessä dummy-demodata
   - Siis, tuotaessa dataa systeemiin, on tuotava myös phenodata ennen 
     normalisointia.
   - Normalisoinnin yhteydessä phenodataan lisätään yksi sarake, jolla kerrotaan
     sirutyyppi, esimerkiksi Affymetrix-sirun versio.
- Laadun tarkistus (QC) sekä Affymetrix että cDNA/Illumina-datalle
- Filtteröintiä parannettu, mm.
   - Ekspressioarvon perusteella filtteröitäessä voidaan määrittää kuinka monella
     sirulla arvon tulee toteutua, jotta arvo läpäisee filtterin.
- Monipuolisemmat tilastolliset työkalut
   - LPE, SAM, etc.
   - Yhdelle sirulle soveltuvat menetelmät (noise-envelope ja Newton)
- Toimiva clusterointi ja useita menetelmiä
   - Etäisyyden voi valita etäisyyksien tai korrelaatioiden joukosta
   - Hierarkkisessa klusteroinnissa tuloksen luotettavuuden voi testata permutaatio-
     testillä
- Pääkomponenttianalyysi
   - Datasetti koodataan uudelleen halutulla määrällä pääkomponentteja
- KNN-menetelmään perustuva "diskriminanttianalyysi"
- Datan kirjoitus ulos GEO:n SOFT-formaatissa
   - Edellyttää, että data on kuvattu hyvin phenodata.txt-tiedostossa, ks.
     http://staff.csc.fi/~jtuimala/SOFT/
- Muut toiminnallisuudet
   - Taulukoiden yhdistäminen geenin nimien perusteella
   - Geenien etsintä korrelaation, nimen, Affymetrix ID:n tai kromosomilokaation 
     perusteella


---------------------------------------------------------------------------------------------------------------------------

R Murskalla

- 64-bittinen versio
   gzip R-2.6.1.tar.gz
   tar xvf R-2.6.1.tar
   cd R-2.6.1
   ./configure
   make
   make check
- Lisäpaketit
  # Basic CRAN packages
  install.packages(c("Matrix", "lme4", "amap", "ape", "flexclust", "kohonen", "e1071", "sma", "fastICA"), repos="http://cran.r-project.org", dependencies = T)
  install.packages(c("XML"), repos="http://cran.r-project.org", dependencies=T)
  install.packages(c("aplpack", "corrgram", "deal", "outliers", "pvclust", "zoo"), repos="http://cran.r-project.org", dependencies=T)
  install.packages(c("pvclust"), repos="http://cran.r-project.org", dependencies=F)

  # Basic Bioconductor packages
  source("http://www.bioconductor.org/biocLite.R")
  biocLite()
  biocLite(c("ctc", "ssize", "LPE", "Ruuid", "graph", "affyQCReport", "GlobalAncova", "impute", "idiogram", "GOstats", "beadarray", "GeneTS", "simpleaffy", "globaltest"))
  biocLite(c("geneplotter", "biomaRt", "lumi", "prada", "siggenes", "plier")) 
  biocLite("cosmo")

  # Annotation packages for Affymetrix
  # Human and yeast
  biocLite(c("hgu133a", "hgu133acdf", "hgu133aprobe", "ygs98", "ygs98cdf", "ygs98probe", "hgu133a2", "hgu133a2cdf", "hgu133a2probe", "hgu133plus2", "hgu133plus2cdf", "hgu133plus2probe"))
  biocLite(c("hs133ahsrefseq", "hs133ahsrefseqcdf", "hs133ahsrefseqprobe", "hs133av2hsrefseq", "hs133av2hsrefseqcdf", "hs133av2hsrefseqprobe", "hs133phsrefseq", "hs133phsrefseqcdf", "hs133phsrefseqprobe"))
  biocLite(c("GO", "KEGG"))
  biocLite("PFAM")

  # Mouse and Rat
  biocLite(c("mgu74a", "mgu74acdf", "mgu74aprobe", "mgu74av2", "mgu74av2cdf", "mgu74av2probe", "moe430a", "moe430acdf", "moe430aprobe",
             "mouse4302", "mouse4302cdf", "mouse4302probe", "mouse430a2", "mouse430a2cdf", "mouse430a2probe", 
             "rae230a", "rae230acdf", "rae230aprobe", "rat2302", "rat2302cdf", "rat2302probe", 
             "rgu34a", "rgu34acdf", "rgu34aprobe"))

  # Other organisms, arabidopsis, C. elegans, drosophila, xenopus
  biocLite(c("ag", "agcdf", "agprobe", "ath1121501", "ath1121501cdf", "ath1121501probe", "celegans", "celeganscdf", "celegansprobe",
             "drosgenome1", "drosgenome1cdf", "drosgenome1probe", "drosophila2", "drosophila2cdf", "drosophila2probe",
             "test3probe", "xenopuslaevis", "xenopuslaeviscdf", "xenopuslaevisprobe", "zebrafish", "zebrafishcdf", "zebrafishprobe",
             "yeast2", "yeast2cdf", "yeast2probe", "test3cdf"))
  biocLite(c("hgu95a", "hgu95acdf", "hgu95aprobe", "hgu95av2", "hgu95av2cdf", "hgu95av2probe"))

  # Alternative annotations
  biocLite(c(
  "hs133ahsensg", "hs133ahsensgcdf", "hs133ahsensgprobe", "hs133av2hsensg", "hs133av2hsensgcdf", "hs133av2hsensgprobe", 
  "hs133phsensg", "hs133phsensgcdf", "hs133phsensgprobe", "mm430mmensg", "mm430mmensgcdf", "mm430mmensgprobe",
  "mm74av1mmensg", "mm74av1mmensgcdf", "mm74av1mmensgprobe", "mm74av2mmensg", "mm74av2mmensgcdf", "mm74av2mmensgprobe",
  "rn230arnensg", "rn230arnensgcdf", "rn230arnensgprobe", "rn230rnensg", "rn230rnensgcdf", "rn230rnensgprobe",
  "rn34arnensg", "rn34arnensgcdf", "rn34arnensgprobe"
  )) 

  # Other organisms with no annotations
  biocLite(c("barley1cdf", "barley1probe", "bovinecdf", "bovineprobe", "canine2cdf", "canine2probe", "caninecdf", "canineprobe",
             "chickencdf", "chickenprobe", "citruscdf", "citrusprobe", "maizecdf", "maizeprobe", "medicagocdf", "medicagoprobe",
             "plasmodiumanophelescdf", "plasmodiumanophelesprobe", "poplarcdf", "poplarprobe", "porcinecdf", "porcineprobe",
             "rhesuscdf", "rhesusprobe", "ricecdf", "riceprobe", "soybeancdf", "soybeanprobe", "sugarcanecdf", "sugarcaneprobe",
             "tomatocdf", "tomatoprobe", "vitisviniferacdf", "vitisviniferaprobe", "wheatcdf", "wheatprobe"))

  # Annotation packages for Agilent
  biocLite(c("hgug4100a", "hgug4101a", "hgug4110b", "hgug4111a", "hgug4112a"))
  biocLite(c("rgug4105a", "rgug4130a"))
  biocLite(c("mgug4104a", "mgug4120a", "mgug4121a", "mgug4122a"))

  # Annotation packages for Illumina
  biocLite(c("illuminaHumanv1", "illuminaHumanv2", "illuminaMousev1", "lumiHumanV1", "lumiHumanV2", "lumiMouseV1", "lumiRatV1"))

  # Extra packages
  # Download: http://addictedtor.free.fr/packages/
  # Needed: fpc, A2R
  # Some annotation packages (Illumina ProbeIDs and Affymetrix exon arrays) are available on a disk.


---------------------------------------------------------------------------------------------------------------------------
