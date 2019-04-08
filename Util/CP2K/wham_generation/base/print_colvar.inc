  &FREE_ENERGY
    &METADYN
      DO_HILLS FALSE
      DELTA_T 300.0
      NT_HILLS 400
      WELL_TEMPERED TRUE
      WW 0.001594 
      &METAVAR
        COLVAR 1
        SCALE 1.0
      &END METAVAR
      &PRINT
        &COLVAR
          &EACH
              MD 1
          &END EACH
          FILENAME =./out.colvar
        &END COLVAR
      &END PRINT
    &END METADYN
  &END FREE_ENERGY
