.section .init
.global main

lui	sp, 0x80005

/* call main */
jal ra, main

/* break */
ebreak //系统自陷

