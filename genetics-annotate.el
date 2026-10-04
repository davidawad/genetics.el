;;; genetics-annotate.el --- Annotation files, risk-allele logic and APOE for genetics.el -*- lexical-binding: t; -*-

;; Copyright (C) 2026 David Awad

;; Author: David Awad <davidawad@protonmail.com>
;; Maintainer: David Awad <davidawad@protonmail.com>
;; Version: 0.1.0
;; Package-Requires: ((emacs "29.1"))
;; Keywords: data, tools
;; URL: https://gitlab.com/davidawad/genetics-el

;; This file is not part of GNU Emacs.

;; Permission is hereby granted, free of charge, to any person obtaining a
;; copy of this software, to deal in it without restriction (MIT licence).

;;; Commentary:

;; Annotations are user-editable JSON or Org files describing SNPs of
;; interest.  JSON: an array of objects with keys rsid, gene, risk_allele,
;; other_allele, effect, magnitude, notes, url, strand, genotypes (a map
;; from genotype to interpretation) and coordinates, a map from build
;; ("GRCh37", "GRCh38") to {"chrom", "pos", "ref"} (1-based, + strand),
;; cited by coordinate_sources.  Org: one heading per SNP with a property
;; drawer (RSID, GENE, RISK_ALLELE, OTHER_ALLELE, EFFECT, MAGNITUDE, URL,
;; STRAND, GT_<genotype> and GRCH37_CHROM/_POS/_REF, GRCH38_CHROM/_POS/_REF
;; properties); the heading body is the notes text.
;;
;; Coordinates let a kit without rsids (a whole-genome VCF) be matched by
;; position on its own build.  When such a kit lists variant sites only, a
;; curated site that is absent is reported as homozygous reference with
;; call source `inferred-ref', always labelled as inferred.  Everything is
;; informational, not medical advice.

;;; Code:

(require 'cl-lib)
(require 'subr-x)
(require 'json)
(require 'genetics-core)

(cl-defstruct (genetics-annotation (:constructor genetics-annotation-create)
                                   (:copier nil))
  "Annotation of one SNP."
  rsid gene risk-allele other-allele effect magnitude notes url strand
  genotypes coordinates sources)

(defun genetics-build-key (build)
  "Return \"37\" or \"38\" for BUILD (e.g. \"GRCh38\", \"hg19\", 37), else nil."
  (let ((b (downcase (format "%s" (or build "")))))
    (cond ((string-match-p "\\`\\(grch\\)?37\\(\\.[0-9p]+\\)?\\'\\|\\`hg19\\'\\|\\`b37\\'" b) "37")
          ((string-match-p "\\`\\(grch\\)?38\\(\\.[0-9p]+\\)?\\'\\|\\`hg38\\'\\|\\`b38\\'" b) "38"))))

(defun genetics--ann-site (file rsid build chrom pos ref)
  "Validate one coordinate of RSID in FILE; return (BUILD-KEY . PLIST).
BUILD, CHROM, POS and REF are the raw values; POS may be a string."
  (let ((key (genetics-build-key build))
        (pos (if (stringp pos) (and (string-match-p "\\`[0-9]+\\'" pos)
                                    (string-to-number pos))
               pos)))
    (unless key
      (genetics--ann-error file "%s: unknown build %S (use GRCh37 or GRCh38)"
                           rsid build))
    (unless (and (stringp chrom) (not (string-empty-p chrom)))
      (genetics--ann-error file "%s: %s coordinate has no chrom" rsid build))
    (unless (and (integerp pos) (> pos 0))
      (genetics--ann-error file "%s: %s coordinate needs a positive pos" rsid
                           build))
    (when (and ref (not (string-match-p "\\`[ACGTacgt]+\\'" ref)))
      (genetics--ann-error file "%s: %s ref %S is not a base" rsid build ref))
    (cons key (list :chrom (genetics-normalize-chrom chrom) :pos pos
                    :ref (and ref (upcase ref))))))

;;;; Loading

(defun genetics--ann-error (file fmt &rest args)
  "Signal an annotation error for FILE with message FMT and ARGS."
  (genetics--error 'genetics-annotation-error "%s: %s" file
                   (apply #'format fmt args)))

(defun genetics--ann-genotypes (alist)
  "Convert ALIST of (GENOTYPE . TEXT) to normalized-key pairs."
  (delq nil
        (mapcar (lambda (c)
                  (let ((k (genetics-genotype-key
                            (upcase (format "%s" (car c))))))
                    (when k (cons k (cdr c)))))
                alist)))

(defun genetics--ann-from-json (obj file index)
  "Build an annotation from decoded JSON OBJ (item INDEX of FILE)."
  (let ((rsid (alist-get 'rsid obj)))
    (unless (and (stringp rsid) (not (string-empty-p rsid)))
      (genetics--ann-error file "entry %d has no rsid" index))
    (genetics-annotation-create
     :rsid (downcase rsid) :gene (alist-get 'gene obj)
     :risk-allele (let ((r (alist-get 'risk_allele obj))) (and r (upcase r)))
     :other-allele (let ((r (alist-get 'other_allele obj))) (and r (upcase r)))
     :effect (alist-get 'effect obj)
     :magnitude (alist-get 'magnitude obj)
     :notes (alist-get 'notes obj) :url (alist-get 'url obj)
     :strand (alist-get 'strand obj)
     :genotypes (genetics--ann-genotypes (alist-get 'genotypes obj))
     :coordinates
     (mapcar (lambda (c)
               (let ((v (cdr c)))
                 (unless (and (consp v) (consp (car v)))
                   (genetics--ann-error file "%s: coordinates.%s must be an object"
                                        rsid (car c)))
                 (genetics--ann-site file rsid (symbol-name (car c))
                                     (let ((ch (alist-get 'chrom v)))
                                       (if (numberp ch) (number-to-string ch) ch))
                                     (alist-get 'pos v) (alist-get 'ref v))))
             (alist-get 'coordinates obj))
     :sources (alist-get 'coordinate_sources obj))))

(defun genetics--load-annotation-json (file)
  "Return the annotations in JSON FILE."
  (let ((data (condition-case err
                  (with-temp-buffer
                    (let ((coding-system-for-read 'utf-8))
                      (insert-file-contents file))
                    (genetics--json-parse (buffer-string)
                                          :object-type 'alist :array-type 'list
                                          :null-object nil :false-object nil))
                (json-error
                 (genetics--ann-error file "invalid JSON (%s)"
                                      (error-message-string err))))))
    (unless (and (listp data)
                 (cl-every (lambda (o) (and (consp o) (consp (car o)))) data))
      (genetics--ann-error file "top level must be an array of objects"))
    (let ((i 0))
      (mapcar (lambda (o) (genetics--ann-from-json o file (cl-incf i))) data))))

(defun genetics--ann-from-org (props body file)
  "Build an annotation from Org PROPS (alist) and BODY text of FILE."
  (let ((rsid (cdr (assoc "RSID" props)))
        (get (lambda (k) (cdr (assoc k props)))))
    (unless rsid
      (genetics--ann-error file "heading without an RSID or rsNNN title"))
    (genetics-annotation-create
     :rsid (downcase rsid) :gene (funcall get "GENE")
     :risk-allele (let ((r (funcall get "RISK_ALLELE"))) (and r (upcase r)))
     :other-allele (let ((r (funcall get "OTHER_ALLELE"))) (and r (upcase r)))
     :effect (funcall get "EFFECT")
     :magnitude (let ((m (funcall get "MAGNITUDE")))
                  (and m (if (string-match-p "\\`[0-9.]+\\'" m)
                             (string-to-number m) m)))
     :notes (let ((n (or (funcall get "NOTES") (string-trim body))))
              (and (not (string-empty-p n)) n))
     :url (funcall get "URL") :strand (funcall get "STRAND")
     :genotypes (genetics--ann-genotypes
                 (cl-loop for (k . v) in props
                          when (string-prefix-p "GT_" k)
                          collect (cons (substring k 3) v)))
     :coordinates
     (cl-loop for b in '("GRCH37" "GRCH38")
              for chrom = (funcall get (concat b "_CHROM"))
              for pos = (funcall get (concat b "_POS"))
              when (or chrom pos)
              collect (genetics--ann-site file rsid b chrom pos
                                          (funcall get (concat b "_REF"))))
     :sources (let ((src (funcall get "COORDINATE_SOURCES")))
                (and src (list src))))))

(defun genetics--load-annotation-org (file)
  "Return the annotations in Org FILE (headings with property drawers)."
  (with-temp-buffer
    (let ((coding-system-for-read 'utf-8))
      (insert-file-contents file))
    (goto-char (point-min))
    (let ((result nil) (props nil) (body nil) (in-drawer nil) (heading nil))
      (cl-flet ((flush ()
                  (when heading
                    (let ((p props))
                      (unless (assoc "RSID" p)
                        (when (string-match "\\<\\(rs[0-9]+\\)\\>" heading)
                          (push (cons "RSID" (match-string 1 heading)) p)))
                      (when (or p (assoc "RSID" p))
                        (push (genetics--ann-from-org
                               p (string-join (nreverse body) "\n") file)
                              result))))
                  (setq props nil body nil in-drawer nil)))
        (while (not (eobp))
          (let ((line (buffer-substring-no-properties
                       (line-beginning-position) (line-end-position))))
            (cond
             ((string-match "\\`\\*+ +\\(.*\\)" line)
              (let ((title (match-string 1 line)))
                (flush)
                (setq heading title)))
             ((string-match-p "\\`[ \t]*:PROPERTIES:[ \t]*\\'" line)
              (setq in-drawer t))
             ((string-match-p "\\`[ \t]*:END:[ \t]*\\'" line)
              (setq in-drawer nil))
             ((and in-drawer
                   (string-match "\\`[ \t]*:\\([A-Za-z0-9_]+\\):[ \t]*\\(.*\\)\\'"
                                 line))
              (push (cons (upcase (match-string 1 line))
                          (string-trim (match-string 2 line)))
                    props))
             (heading (push line body))))
          (forward-line 1))
        (flush))
      (nreverse result))))

(defun genetics-load-annotation-file (file)
  "Return the list of annotations read from FILE (.json or .org)."
  (unless (file-readable-p file)
    (genetics--ann-error file "file not found or unreadable"))
  (pcase (downcase (or (file-name-extension file) ""))
    ("json" (genetics--load-annotation-json file))
    ("org" (genetics--load-annotation-org file))
    (ext (genetics--ann-error file "unsupported extension %S (use .json or .org)"
                              ext))))

(defvar genetics--annotation-cache nil
  "Cons of (SIGNATURE . HASH-TABLE) for the loaded annotations.")

(defvar genetics--annotation-position-cache nil
  "Cons of (HASH-TABLE . POSITION-TABLE) indexing annotations by site.
POSITION-TABLE maps \"BUILD:CHROM:POS\" to an annotation.")

(defun genetics--annotation-signature ()
  "Return a value that differs once annotation files are modified."
  (mapcar (lambda (f)
            (list f (file-attribute-modification-time (file-attributes f))
                  (file-attribute-size (file-attributes f))))
          genetics-annotation-files))

(defun genetics-annotations ()
  "Return a hash table of rsid to annotation for `genetics-annotation-files'."
  (let ((sig (genetics--annotation-signature)))
    (unless (and genetics--annotation-cache
                 (equal sig (car genetics--annotation-cache)))
      (let ((table (make-hash-table :test 'equal)))
        (dolist (f genetics-annotation-files)
          (dolist (a (genetics-load-annotation-file f))
            (puthash (genetics-annotation-rsid a) a table)))
        (setq genetics--annotation-cache (cons sig table))))
    (cdr genetics--annotation-cache)))

(defun genetics-reload-annotations ()
  "Discard the annotation cache so files are re-read on next use."
  (interactive)
  (setq genetics--annotation-cache nil)
  (hash-table-count (genetics-annotations)))

(defun genetics-annotation (rsid)
  "Return the annotation for RSID, or nil."
  (gethash rsid (genetics-annotations)))

(defun genetics-annotation-site (ann build)
  "Return the plist (:chrom :pos :ref) of ANN on BUILD, or nil."
  (cdr (assoc (genetics-build-key build) (genetics-annotation-coordinates ann))))

(defun genetics-annotation-at (build chrom pos)
  "Return the annotation whose BUILD coordinate is CHROM:POS, or nil."
  (let ((table (genetics-annotations)))
    (unless (eq (car genetics--annotation-position-cache) table)
      (let ((pt (make-hash-table :test 'equal)))
        (maphash (lambda (_r ann)
                   (dolist (c (genetics-annotation-coordinates ann))
                     (puthash (format "%s:%s:%d" (car c)
                                      (plist-get (cdr c) :chrom)
                                      (plist-get (cdr c) :pos))
                              ann pt)))
                 table)
        (setq genetics--annotation-position-cache (cons table pt))))
    (let ((key (genetics-build-key build)))
      (when key
        (gethash (format "%s:%s:%d" key chrom pos)
                 (cdr genetics--annotation-position-cache))))))

(defun genetics-snp-annotation (snp &optional build annotations)
  "Return the annotation of SNP by rsid, or by its position on BUILD.
ANNOTATIONS is the rsid table (default `genetics-annotations')."
  (or (gethash (genetics-snp-rsid snp) (or annotations (genetics-annotations)))
      (and build (genetics-annotation-at build (genetics-snp-chrom snp)
                                         (genetics-snp-pos snp)))))

;;;; Resolving curated sites in a kit

(defconst genetics--builtin-sites
  '(("rs429358" ("37" :chrom "19" :pos 45411941 :ref "T")
                ("38" :chrom "19" :pos 44908684 :ref "T"))
    ("rs7412" ("37" :chrom "19" :pos 45412079 :ref "C")
              ("38" :chrom "19" :pos 44908822 :ref "C")))
  "Coordinates of the APOE SNPs, used when no annotation file has them.
Verified 2026-10-03 against NCBI dbSNP and Ensembl (see
annotations/genetics-curated.json).")

(defun genetics-rsid-site (rsid build)
  "Return the (:chrom :pos :ref) plist of RSID on BUILD, or nil.
Coordinates come from the annotations, then `genetics--builtin-sites'."
  (let ((ann (genetics-annotation rsid)))
    (or (and ann (genetics-annotation-site ann build))
        (cdr (assoc (genetics-build-key build)
                    (cdr (assoc rsid genetics--builtin-sites)))))))

(defun genetics--autosome-p (chrom)
  "Return non-nil if CHROM is one of chromosomes 1-22."
  (and (string-match-p "\\`[0-9]+\\'" chrom)
       (<= 1 (string-to-number chrom) 22)))

(defun genetics-kit-infer-ref (kit rsid site)
  "Return an inferred homozygous-reference call for RSID at SITE in KIT.
Return nil unless KIT is a variant-only WGS kit (`absent-means-ref'),
SITE has a reference base on an autosome the kit has records for, and no
non-reference deletion in KIT spans the site."
  (let ((chrom (plist-get site :chrom)) (pos (plist-get site :pos))
        (ref (plist-get site :ref)))
    (when (and (eq (genetics-kit-ref-calls kit) 'absent-means-ref)
               ref (genetics--autosome-p chrom)
               (member chrom (genetics-kit-chromosomes kit))
               (not (genetics-kit-site-spanned-p kit chrom pos)))
      (genetics-snp-create rsid chrom pos (concat ref ref) ref nil
                           'inferred-ref))))

(defun genetics-kit-resolve (kit rsid)
  "Return the call for RSID in KIT as a `genetics-snp', or nil.
Looks up RSID itself; then, using the curated coordinates on the kit's
build, the record at that position; then, for variant-only WGS kits, an
inferred homozygous-reference call (call source `inferred-ref').  Kits
served by genome-cli resolve all of this in genome-cli."
  (if (genetics-kit-backend kit)
      (genetics-kit-get kit rsid)
    (or (genetics-kit-get kit rsid)
        (let ((site (genetics-rsid-site rsid (genetics-kit-build kit))))
          (when site
            (or (genetics-kit-at kit (plist-get site :chrom)
                                 (plist-get site :pos) (plist-get site :ref))
                (genetics-kit-infer-ref kit rsid site)))))))

(defun genetics-annotated-snps (kit)
  "Return a list of (ANNOTATION . SNP) for annotated rsids present in KIT.
Sites are resolved with `genetics-kit-resolve', so a kit without rsids is
matched by position and SNP may be an inferred reference call."
  (let (hits)
    (maphash (lambda (rsid ann)
               (let ((snp (genetics-kit-resolve kit rsid)))
                 (when snp (push (cons ann snp) hits))))
             (genetics-annotations))
    (sort hits (lambda (a b)
                 (let ((ga (or (genetics-annotation-gene (car a)) ""))
                       (gb (or (genetics-annotation-gene (car b)) "")))
                   (if (string= ga gb)
                       (string< (genetics-annotation-rsid (car a))
                                (genetics-annotation-rsid (car b)))
                     (string< ga gb)))))))

;;;; Risk allele logic

(defun genetics-risk-assess (genotype risk &optional other)
  "Assess GENOTYPE against RISK allele (OTHER is the non-risk allele).
Return a plist with :copies (nil if no-call or unusable), :flag and
:flipped-copies.  :flag is nil, `no-call', `ambiguous' (palindromic SNP),
`strand-flip' (only complementary alleles seen) or `unexpected'.  The
genotype is never silently flipped."
  (let* ((alleles (genetics-alleles genotype))
         (risk (and risk (upcase risk)))
         (other (and other (upcase other)))
         (crisk (genetics-complement risk))
         (cother (and other (genetics-complement other)))
         (copies (cl-count risk alleles :test #'equal)))
    (cond
     ((null alleles) (list :copies nil :flag 'no-call))
     ((null risk) (list :copies nil :flag nil))
     ;; palindromic: the SNP's own alleles are complements of each other
     ((or (and other crisk (equal other crisk))
          (and (null other) crisk (member risk alleles) (member crisk alleles)))
      (list :copies copies :flag 'ambiguous))
     ;; none of the observed alleles is risk/other, but complements are
     ((and (zerop copies) (not (member other alleles))
           (cl-some (lambda (a) (or (equal a crisk) (equal a cother))) alleles)
           (cl-every (lambda (a) (or (equal a crisk) (equal a cother))) alleles))
      (list :copies 0 :flag 'strand-flip
            :flipped-copies (cl-count crisk alleles :test #'equal)))
     ((and other
           (cl-notevery (lambda (a) (or (equal a risk) (equal a other)))
                        alleles))
      (list :copies copies :flag 'unexpected))
     (t (list :copies copies :flag nil)))))

(defun genetics-flag-text (flag flipped)
  "Return text explaining assessment FLAG (FLIPPED is the flipped copy count)."
  (pcase flag
    ('no-call "No call at this position.")
    ('ambiguous "Palindromic SNP (A/T or C/G): the strand cannot be verified from the genotype, treat the copy count with caution.")
    ('strand-flip (format "Possible strand flip: the genotype only contains complements of the expected alleles. Not flipped automatically; if the data are on the opposite strand the risk allele copies would be %s." flipped))
    ('unexpected "Genotype contains alleles other than the expected risk/other alleles; check strand and build.")
    (_ nil)))

(defun genetics--copies-text (copies risk zyg)
  "Return generic text for COPIES of RISK allele given zygosity ZYG."
  (pcase copies
    (0 (format "No copies of the %s allele." risk))
    (1 (format "One copy of the %s allele (%s)." risk
               (if (eq zyg 'hemizygous) "hemizygous" "heterozygous")))
    (2 (format "Two copies of the %s allele (homozygous)." risk))
    (_ (format "%s copies of the %s allele." copies risk))))

(defun genetics-assess (ann genotype)
  "Assess GENOTYPE against annotation ANN; return a plist.
Keys: :genotype :copies :flag :flipped-copies :interpretation :strand."
  (let* ((a (genetics-risk-assess genotype (genetics-annotation-risk-allele ann)
                                  (genetics-annotation-other-allele ann)))
         (flag (plist-get a :flag))
         (copies (plist-get a :copies))
         (custom (cdr (assoc (genetics-genotype-key genotype)
                             (genetics-annotation-genotypes ann))))
         (risk (genetics-annotation-risk-allele ann))
         (text
          (cond
           ((eq flag 'no-call) "No call; nothing can be said about this SNP.")
           ((eq flag 'strand-flip)
            (genetics-flag-text flag (plist-get a :flipped-copies)))
           ((eq flag 'unexpected) (genetics-flag-text flag nil))
           (t (string-join
               (delq nil
                     (list (or custom
                               (and copies risk
                                    (genetics--copies-text
                                     copies risk (genetics-zygosity genotype)))
                               (genetics-annotation-effect ann))
                           (genetics-flag-text flag nil)))
               " ")))))
    (append (list :genotype genotype :interpretation text
                  :strand (genetics-annotation-strand ann))
            a)))

(defun genetics-assess-snp (ann snp)
  "Assess call SNP against annotation ANN, as `genetics-assess'.
An inferred call has `genetics-inferred-ref-text' prepended to its
interpretation and :call-source `inferred-ref'."
  (let ((a (genetics-assess ann (genetics-snp-genotype snp))))
    (if (genetics-snp-inferred-p snp)
        (plist-put (plist-put a :interpretation
                              (concat genetics-inferred-ref-text " "
                                      (plist-get a :interpretation)))
                   :call-source 'inferred-ref)
      a)))

;;;; APOE

(defconst genetics-apoe-snps '("rs429358" "rs7412")
  "The two SNPs that define the APOE epsilon haplotypes.")

(defun genetics--allele-count (genotype allele)
  "Count ALLELE in GENOTYPE."
  (cl-count allele (genetics-alleles genotype) :test #'equal))

(defun genetics-apoe-interpret (g429358 g7412)
  "Return the APOE result for genotypes G429358 (rs429358) and G7412 (rs7412).
Assumes + strand (rs429358 T/C, rs7412 C/T): e2 = T+T, e3 = T+C,
e4 = C+C, e1 = C+T.  Returns a plist with :diplotype, :status (`ok',
`ambiguous', `unusual', `strand-flip', `no-call' or `incomplete') and
:description."
  (let ((a1 (genetics-alleles g429358)) (a2 (genetics-alleles g7412)))
    (cond
     ((or (null a1) (null a2))
      (list :diplotype nil :status 'no-call
            :description "APOE cannot be determined: a no-call at rs429358 or rs7412."))
     ((or (/= (length a1) 2) (/= (length a2) 2))
      (list :diplotype nil :status 'incomplete
            :description "APOE needs diploid genotypes at both SNPs."))
     ((or (cl-notevery (lambda (x) (member x '("C" "T"))) a1)
          (cl-notevery (lambda (x) (member x '("C" "T"))) a2))
      (if (and (cl-every (lambda (x) (member x '("A" "G"))) a1)
               (cl-every (lambda (x) (member x '("A" "G"))) a2))
          (list :diplotype nil :status 'strand-flip
                :description "Possible strand flip: genotypes use A/G instead of T/C. Not flipped automatically; verify the strand before interpreting APOE.")
        (list :diplotype nil :status 'incomplete
              :description "Unexpected alleles at rs429358/rs7412; APOE not determined.")))
     (t
      (let* ((c1 (genetics--allele-count g429358 "C"))
             (t2 (genetics--allele-count g7412 "T"))
             (res (pcase (cons c1 t2)
                    ('(0 . 0) '("e3/e3" ok))
                    ('(0 . 1) '("e2/e3" ok))
                    ('(0 . 2) '("e2/e2" ok))
                    ('(1 . 0) '("e3/e4" ok))
                    ('(1 . 1) '("e2/e4" ambiguous))
                    ('(1 . 2) '("e1/e2" unusual))
                    ('(2 . 0) '("e4/e4" ok))
                    ('(2 . 1) '("e1/e4" unusual))
                    ('(2 . 2) '("e1/e1" unusual))))
             (dip (car res)) (status (cadr res)))
        (list :diplotype dip :status status
              :description
              (pcase status
                ('ambiguous "e2/e4 (most likely) or e1/e3 (very rare): the phase of a double heterozygote cannot be resolved from unphased genotypes, so e2/e4 is assumed.")
                ('unusual (format "%s is a rare combination; confirm with a clinical-grade assay (e1 is very rare)." dip))
                (_ (genetics--apoe-blurb dip)))))))))

(defun genetics--apoe-blurb (dip)
  "Return a short informational statement for APOE diplotype DIP."
  (format "%s. %s" dip
          (pcase dip
            ("e3/e3" "The most common genotype; baseline late-onset Alzheimer's disease risk.")
            ("e2/e3" "One e2 allele; associated with somewhat lower late-onset Alzheimer's risk than e3/e3.")
            ("e2/e2" "Two e2 alleles; associated with lower Alzheimer's risk and with type III hyperlipoproteinemia in some carriers.")
            ("e3/e4" "One e4 allele; associated with higher late-onset Alzheimer's risk than e3/e3 (not deterministic).")
            ("e4/e4" "Two e4 alleles; associated with substantially higher late-onset Alzheimer's risk (not deterministic).")
            (_ ""))))

(defun genetics-apoe-for-kit (kit)
  "Return the APOE result plist for KIT, or nil if the SNPs are absent.
The SNPs are found with `genetics-kit-resolve'.  The plist gains
:inferred, the list of APOE rsids whose call was inferred rather than
observed, and the description says so."
  (let ((a (genetics-kit-resolve kit "rs429358"))
        (b (genetics-kit-resolve kit "rs7412")))
    (when (and a b)
      (let ((res (genetics-apoe-interpret (genetics-snp-genotype a)
                                          (genetics-snp-genotype b)))
            (inferred (delq nil (list (and (genetics-snp-inferred-p a) "rs429358")
                                      (and (genetics-snp-inferred-p b) "rs7412")))))
        (if inferred
            (plist-put (plist-put res :inferred inferred) :description
                       (format "%s %s inferred homozygous reference (absent from a variant-only WGS VCF), not observed."
                               (plist-get res :description)
                               (string-join inferred " and ")))
          res)))))

(provide 'genetics-annotate)
;;; genetics-annotate.el ends here
