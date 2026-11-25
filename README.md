# Classification des SNV constitutionnels — ACGS 2026

Application d’aide à l’interprétation des substitutions codantes ou d’épissage
sur les chromosomes nucléaires primaires GRCh38. Le navigateur local et la
commande VCF utilisent le même moteur Python. Les fichiers sont lus sur le poste
de travail. Aucune connexion Internet ni installation de dépendances ne sont
nécessaires à l’exécution.

## Démarrer

**Windows x64.** Extraire l’archive dans un dossier utilisateur, puis ouvrir
`launch_windows.cmd`. L’archive inclut CPython x64 avec sa licence. Le navigateur
ouvre `http://127.0.0.1:8765`. Laisser la console ouverte pendant l’utilisation.
Les champs navigateur restent dans la session. L’export JSON s’effectue à la
demande. Utiliser le bouton de validation après avoir relu les preuves.

Pour essayer l’interface, saisir `chr1:123:A:G` avec la transmission AD et coller
`examples/annotation_a_coller.json` dans les annotations. Le résultat simulé est
VUS à 3 points (PM2 et PP3). Ajouter le fichier `examples/ps3.tsv` donne LP à
7 points si le contexte est inchangé. Les exemples sont entièrement simulés.

**Linux.** Depuis le dossier extrait, utiliser Python 3.10 ou plus récent.

```sh
python3 -m acmg_snv --help
python3 -m acmg_snv web --no-browser
python3 -m unittest discover -s tests -v
```

La commande `web` est également disponible sous Linux. Le serveur Web accepte
uniquement les requêtes locales et n’est pas destiné à un déploiement partagé.

## Configurer les données du laboratoire

Seuls le catalogue ClinGen et les jeux de démonstration sont inclus. Préparer les
autres ressources dans un dossier accessible en lecture. Utiliser la même
version de GRCh38 pour la référence, l’annotation et les variants.

```sh
python3 -m acmg_snv index-fasta reference.fa
python3 -m acmg_snv index-mane MANE.GRCh38.refseq_genomic.gff.gz mane.sqlite
python3 -m acmg_snv index annotations.jsonl annotations.sqlite
```

La référence est un FASTA **non compressé** avec index `.fai`. Réutiliser un index
existant ou le créer avec `index-fasta`. L’index MANE utilise le **GFF3 RefSeq
génomique GRCh38** de la version MANE du laboratoire. Le couple FASTA/index
MANE permet de vérifier REF et de calculer la conséquence et la nomenclature
pour les transcrits MANE couverts. Un annotateur externe peut fournir des
conséquences ou des éléments d’épissage supplémentaires.

L’index d’annotations est facultatif. Il accepte du JSONL, JSON ou TSV préparé à
partir des annotations du pipeline, avec une ligne par allèle et transcrit.
Les quatre champs `chrom`, `pos`, `ref`, `alt` identifient l’allèle
(orientation génomique, position à partir de 1). Ajouter `gene`, `transcript`
(avec version), `mane`, `consequence` (terme Sequence Ontology), `source` et
`version` quand ils sont disponibles. Les champs de fréquences peuvent aussi
être fournis sur une ligne par allèle, sans gène ni transcrit.

Dans le navigateur ou via `classify --annotation`, coller également un objet
JSON, un TSV normalisé, ou l’en-tête `##INFO=<ID=CSQ,...>` suivi de la ligne VCF
annotée avec VEP. Les colonnes et les allèles VEP sont interprétés d’après
l’en-tête. Les formats propriétaires et les bases brutes doivent être
convertis. La référence FASTA brute et les bases gnomAD et prédicteurs peuvent
provenir de votre pipeline local. Une annotation ClinVar sert de piste de
revue et n’attribue aucun point à la seule assertion publiée.

Démarrer Windows depuis l’invite de commandes avec les ressources configurées.

```bat
launch_windows.cmd --reference reference.fa --mane-db mane.sqlite ^
  --annotation-db annotations.sqlite
```

Sur Linux, passer les mêmes options à `python3 -m acmg_snv web`, `classify`
ou `vcf`. Depuis le dossier racine, le lanceur Windows accepte également
les sous-commandes Python.

```bat
launch_windows.cmd index-fasta reference.fa
```

## Traitement VCF et SLURM

```sh
python3 -m acmg_snv index examples/annotations.jsonl demo.sqlite
python3 -m acmg_snv vcf examples/variants.vcf -o demo-result.vcf.gz \
  --context examples/context.json --annotation-db demo.sqlite
```

`vcf` accepte le VCF en texte, gzip ou BGZF. Il conserve les colonnes et
ajoute par ALT `ACMG_CLASS`, `ACMG_SCORE` et `ACMG_STATUS`. Les ALT hors
périmètre ou sans données exploitées sont identifiables. Les sorties `.gz`
sont en BGZF. Le JSONL (par défaut `SORTIE.jsonl`) contient les justifications,
les données manquantes, la force, les points et les alertes. Un flux explicite
sur entrée ou sortie standard utilise `-`. Avec la sortie standard, fournir
`--report` vers un autre fichier. Les modes de transmission AD, AR, XL, XLD,
XLR, YL et OTHER sont acceptés. La classe est suspendue si la transmission
manque. Pour des maladies différentes selon l’allèle, fournir `--contexts`
avec un objet JSON dont les clés sont les coordonnées normalisées.

```sh
sbatch slurm.sh variants.vcf.gz -o resultat.vcf.gz \
  --contexts contextes_par_allele.json --annotation-db annotations.sqlite \
  --reference reference.fa --mane-db mane.sqlite
```

La commande hors ligne `vcep --gene BRCA1` expose le détail du catalogue.

Le script SLURM demande un seul CPU et réutilise
les options de la commande `vcf`. Exécuter depuis la racine de l’application.
Adapter les ressources de la tâche à la taille du VCF et au stockage du cluster.

## Preuves et interprétation

- Les forces pathogènes supporting, moderate, strong et very strong valent
  +1, +2, +4 et +8. Les niveaux bénins correspondants ont le signe opposé.
  Le score SNV de 2026 donne P à partir de 10, LP de 6 à 9, VUS de 0 à 5,
  LB de −1 à −6 et B à partir de −7 vers le bas. Il faut des preuves
  corroborantes pour classer hors VUS. BA1 est autonome et n’a pas de valeur
  numérique. BS1 fort peut être autonome pour une affection AD hautement
  pénétrante dans le référentiel général.
- Les critères calculés sont des propositions à partir des données explicites.
  Les preuves de littérature et les essais ne sont pas interprétés à partir
  de texte libre. Une même preuve ne peut être comptée deux fois.
  Les plafonds et incompatibilités documentés sont indiqués dans le rapport.
  Un désaccord entre preuves pathogènes et bénignes demande une revue motivée.
- Les propositions générales couvrent notamment les fréquences
  BA1/BS1/PM2, les prédictions PP3/BP4/BP7, des branches de PVS1, les
  comparaisons PS1/PM5 documentées et les données cliniques structurées
  PS2/PM6, PM3, PP1 et PP4. Autres situations à documenter en revue.
- Le résumé clinique et la ségrégation en texte libre sont conservés dans la
  fiche. Pour une proposition cotée, ajouter des observations structurées,
  ou sélectionner le critère et sa force avec une justification et une
  source. Utiliser des identifiants de preuve stables pour repérer les
  observations réutilisées.
- Le fichier fonctionnel PS3 est un TSV ou CSV UTF-8 comportant
  `variant`, `strength` et `reason`. `source` et `evidence_id` sont facultatifs.
  Une ligne unique par allèle, issue d’un essai protéique spécifique validé,
  a priorité sur les preuves PS3 importées de la littérature. Les essais
  d’épissage ARN relèvent de l’évaluation PVS1(RNA) ou BP7(RNA) selon S4.
- Pour un gène associé à plusieurs maladies, renseigner la maladie et le
  mécanisme pertinents. Quand plusieurs gènes ou transcrits sont candidats,
  sélectionner le contexte clinique. MANE Select est prioritaire par défaut,
  avec MANE Plus Clinical et autres transcrits disponibles selon les données.
- `provisional` signifie proposition non validée. Les cas `needs_context`,
  `not_evaluable`, `out_of_scope` et `review_required` peuvent avoir une classe
  vide. `approved` correspond à la validation explicitement saisie par
  l’utilisateur dans cette session.

## ClinGen VCEP

L’instantané du CSpec contient les documents détaillés des 127 entrées
`Released` récupérées le 24 novembre 2026. Les versions explicitement
remplacées sont exclues de la sélection active. Les spécifications sont
identifiées par gène, maladie et transmission. Utiliser si possible
l’identifiant MONDO du contexte et confirmer la spécification dans la
fiche. L’interface affiche les règles CSpec disponibles hors ligne et
contrôle les forces explicitement exclues.

Lorsqu’un VCEP s’applique, ses règles priment sur les propositions génériques.
La majorité de ces règles requiert une lecture et une cotation par le
biologiste. Le profil MYH7 GN002 v2.0 propose BA1 et BS1 dans les contextes
MONDO et AD couverts, lorsque les données FAF gnomAD qualifiées et la
confirmation de comparabilité des cohortes sont explicites. La décision
VCEP reste soumise à la validation du profil. Une version de catalogue ou
des profils de seuils partiels validés localement sont chargeables avec
`--vcep-dir`.

## Références

Durkie et al., *Association for Clinical Genomic Science Best Practice
Guidelines for Variant Classification in Rare Disease*, 2026,
recommandation 6, figures 1 et 2, suppléments S2, S3 et S4.
[DOI 10.1136/jmg-2026-111835](https://doi.org/10.1136/jmg-2026-111835)

[ClinGen CSpec](https://cspec.genome.network/cspec/ui/svi/srvc) ·
[Catalogue ClinGen](https://cspec.genome.network/cspec/api/svis) ·
[Conditions d’utilisation des données ClinGen](https://clinicalgenome.org/about/terms-of-use/)

[MANE GRCh38](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/) ·
[Python 3.14.10 Windows embarqué](https://www.python.org/downloads/release/python-31410/)
