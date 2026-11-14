/* Synthetic SEER layout for PyHealth tests. Positions are invented. */
data seer;
  infile 'incidence/*.TXT' lrecl=80 pad;
  input
    @ 1    PUBCSNUM           $char8.  /* Patient ID */
    @ 9    REG                $char10.  /* SEER registry */
    @ 19   MAR_STAT           $char1.  /* Marital status at diagnosis */
    @ 20   RACE1V             $char2.  /* Race recode (W, B, AI, API) */
    @ 22   SEX                $char1.  /* Sex */
    @ 23   AGE_DX             $char3.  /* Age at diagnosis */
    @ 26   SEQ_NUM            $char2.  /* Sequence number */
    @ 28   YEAR_DX            $char4.  /* Year of diagnosis */
    @ 32   PRIMSITE           $char3.  /* Primary Site */
    @ 35   LATERAL            $char1.  /* Laterality */
    @ 36   HISTO3V            $char4.  /* Histologic Type ICD-O-3 */
    @ 40   BEHO3V             $char1.  /* Behavior code ICD-O-3 */
    @ 41   GRADE              $char1.  /* Grade */
    @ 42   EOD10_PN           $char2.  /* Regional nodes positive (1988+) */
    @ 44   EOD10_NE           $char2.  /* Regional nodes examined (1988+) */
    @ 46   CSTUMSIZ           $char3.  /* CS tumor size (2004+) */
    @ 49   CSEXTEN            $char3.  /* CS extension (2004+) */
    @ 52   CSLYMPHN           $char3.  /* CS lymph nodes (2004+) */
    @ 55   CSMETSDX           $char2.  /* CS mets at dx (2004+) */
    @ 57   RX_SUMM_SURG       $char2.  /* RX Summ--surg prim site */
    @ 59   REC_NO             $char2.  /* Record number */
    @ 61   FIRSTPRM           $char1.  /* First malignant primary indicator */
    @ 62   ST_CNTY            $char5.  /* State-county */
    @ 67   STAT_REC           $char1.  /* Vital status recode (study cutoff used) */
    @ 68   VSRTSADX           $char1.  /* SEER cause of death classification */
    @ 69   ODTHCLASS          $char1.  /* SEER other cause of death classification */
    @ 70   SRV_TIME_MON       $char4.  /* Survival months */
    @ 74   SRV_TIME_MON_FLAG  $char1.  /* Survival months flag */
    @ 75   DAJCCSTG           $char2.  /* Derived AJCC Stage Group, 6th ed (2004+) */
    @ 77   SITERWHO           $char5.  /* Site recode ICD-O-3/WHO 2008 */
  ;
run;
