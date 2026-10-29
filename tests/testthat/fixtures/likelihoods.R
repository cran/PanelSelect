# Frozen independent reference values; do not regenerate from production code.
# Decimal literals retain 17 significant digits for exact numeric round trips.
list(fixed = list(`probit_re -1` = list(theta = c(`outcome.(Intercept)` = 0.25, 
outcome.x = -0.20000000000000001, sd_outcome = 0.40000000000000002
), panel = list(x = structure(c(1, 1, 1, -0.40000000000000002, 
0.29999999999999999, 0.69999999999999996), dim = 3:2, dimnames = list(
    NULL, c("(Intercept)", "x"))), w = structure(c(1, 1, 1, 0.5, 
-0.29999999999999999, 0.10000000000000001), dim = 3:2, dimnames = list(
    NULL, c("(Intercept)", "w"))), y = c(1, 0, 0), d = c(1, 1, 
1)), outer = 7L, inner = 13L, loglik = -2.147387547785971, score = c(-0.92511074769782331, 
-1.0828359105695453, -0.2283496767056343)), `probit_re 1` = list(
    theta = c(`outcome.(Intercept)` = 0.25, outcome.x = -0.20000000000000001, 
    sd_outcome = 0.40000000000000002), panel = list(x = structure(c(1, 
    1, 1, -0.40000000000000002, 0.29999999999999999, 0.69999999999999996
    ), dim = 3:2, dimnames = list(NULL, c("(Intercept)", "x"))), 
        w = structure(c(1, 1, 1, 0.5, -0.29999999999999999, 0.10000000000000001
        ), dim = 3:2, dimnames = list(NULL, c("(Intercept)", 
        "w"))), y = c(1, 0, 0), d = c(1, 1, 1)), outer = 7L, 
    inner = 13L, loglik = -2.147387547785971, score = c(-0.92511074769782331, 
    -1.0828359105695453, -0.2283496767056343)), `poisson_re -1` = list(
    theta = c(`outcome.(Intercept)` = 0.25, outcome.x = -0.20000000000000001, 
    sd_outcome = 0.40000000000000002), panel = list(x = structure(c(1, 
    1, 1, -0.40000000000000002, 0.29999999999999999, 0.69999999999999996
    ), dim = 3:2, dimnames = list(NULL, c("(Intercept)", "x"))), 
        w = structure(c(1, 1, 1, 0.5, -0.29999999999999999, 0.10000000000000001
        ), dim = 3:2, dimnames = list(NULL, c("(Intercept)", 
        "w"))), y = c(1, 0, 2), d = c(1, 1, 1)), outer = 7L, 
    inner = 13L, loglik = -4.0555286191414437, score = c(-0.56878630094015104, 
    0.43557871346055077, -0.76813534157896979)), `poisson_re 1` = list(
    theta = c(`outcome.(Intercept)` = 0.25, outcome.x = -0.20000000000000001, 
    sd_outcome = 0.40000000000000002), panel = list(x = structure(c(1, 
    1, 1, -0.40000000000000002, 0.29999999999999999, 0.69999999999999996
    ), dim = 3:2, dimnames = list(NULL, c("(Intercept)", "x"))), 
        w = structure(c(1, 1, 1, 0.5, -0.29999999999999999, 0.10000000000000001
        ), dim = 3:2, dimnames = list(NULL, c("(Intercept)", 
        "w"))), y = c(1, 0, 2), d = c(1, 1, 1)), outer = 7L, 
    inner = 13L, loglik = -4.0555286191414437, score = c(-0.56878630094015104, 
    0.43557871346055077, -0.76813534157896979)), `pln_re -1` = list(
    theta = c(`outcome.(Intercept)` = 0.25, outcome.x = -0.20000000000000001, 
    sd_outcome = 0.40000000000000002, sd_error = 0.5), panel = list(
        x = structure(c(1, 1, 1, -0.40000000000000002, 0.29999999999999999, 
        0.69999999999999996), dim = 3:2, dimnames = list(NULL, 
            c("(Intercept)", "x"))), w = structure(c(1, 1, 1, 
        0.5, -0.29999999999999999, 0.10000000000000001), dim = 3:2, dimnames = list(
            NULL, c("(Intercept)", "w"))), y = c(1, 0, 2), d = c(1, 
        1, 1)), outer = 7L, inner = 13L, loglik = -4.1679740304524877, 
    score = c(-0.63647982123000002, 0.29177786205773393, -0.5900069360714344, 
    -0.47032686494040721)), `pln_re 1` = list(theta = c(`outcome.(Intercept)` = 0.25, 
outcome.x = -0.20000000000000001, sd_outcome = 0.40000000000000002, 
sd_error = 0.5), panel = list(x = structure(c(1, 1, 1, -0.40000000000000002, 
0.29999999999999999, 0.69999999999999996), dim = 3:2, dimnames = list(
    NULL, c("(Intercept)", "x"))), w = structure(c(1, 1, 1, 0.5, 
-0.29999999999999999, 0.10000000000000001), dim = 3:2, dimnames = list(
    NULL, c("(Intercept)", "w"))), y = c(1, 0, 2), d = c(1, 1, 
1)), outer = 7L, inner = 13L, loglik = -4.1679740304524877, score = c(-0.63647982123000002, 
0.29177786205773393, -0.5900069360714344, -0.47032686494040721
)), `selection_gaussian -1` = list(theta = c(`selection.(Intercept)` = -0.14999999999999999, 
selection.w = 0.29999999999999999, `outcome.(Intercept)` = 0.25, 
outcome.x = -0.20000000000000001, sd_outcome = 0.40000000000000002, 
sd_selection = 0.45000000000000001, cor_id = 0.34999999999999998, 
sd_error = 0.5, cor_error = -0.25), panel = list(x = structure(c(1, 
1, 1, -0.40000000000000002, 0.29999999999999999, 0.69999999999999996
), dim = 3:2, dimnames = list(NULL, c("(Intercept)", "x"))), 
    w = structure(c(1, 1, 1, 0.5, -0.29999999999999999, 0.10000000000000001
    ), dim = 3:2, dimnames = list(NULL, c("(Intercept)", "w"))), 
    y = c(0.20000000000000001, NA, -0.40000000000000002), d = c(1, 
    0, 1)), outer = 7L, inner = 13L, loglik = -3.2030583566145157, 
    score = c(0.7767826046103794, 0.67663369339234702, -0.96907213930647795, 
    -0.99355232129138704, -1.1513826722396787, -0.41808273699294141, 
    -0.11464816893670188, -1.8726644745970971, -0.51043356201753942
    )), `selection_gaussian 1` = list(theta = c(`selection.(Intercept)` = -0.14999999999999999, 
selection.w = 0.29999999999999999, `outcome.(Intercept)` = 0.25, 
outcome.x = -0.20000000000000001, sd_outcome = 0.40000000000000002, 
sd_selection = 0.45000000000000001, cor_id = -0.34999999999999998, 
sd_error = 0.5, cor_error = 0.25), panel = list(x = structure(c(1, 
1, 1, -0.40000000000000002, 0.29999999999999999, 0.69999999999999996
), dim = 3:2, dimnames = list(NULL, c("(Intercept)", "x"))), 
    w = structure(c(1, 1, 1, 0.5, -0.29999999999999999, 0.10000000000000001
    ), dim = 3:2, dimnames = list(NULL, c("(Intercept)", "w"))), 
    y = c(0.20000000000000001, NA, -0.40000000000000002), d = c(1, 
    0, 1)), outer = 7L, inner = 13L, loglik = -3.3242023130920253, 
    score = c(0.71596462872460975, 0.65801702455273103, -1.2986532265567394, 
    -1.084984234056946, -0.59056464292347211, -0.23119515326008361, 
    -0.19030798482973704, -2.0813724478868987, -0.40910754752373635
    )), `selection_probit -1` = list(theta = c(`selection.(Intercept)` = -0.14999999999999999, 
selection.w = 0.29999999999999999, `outcome.(Intercept)` = 0.25, 
outcome.x = -0.20000000000000001, sd_outcome = 0.40000000000000002, 
sd_selection = 0.45000000000000001, cor_id = 0.34999999999999998, 
cor_error = -0.25), panel = list(x = structure(c(1, 1, 1, -0.40000000000000002, 
0.29999999999999999, 0.69999999999999996), dim = 3:2, dimnames = list(
    NULL, c("(Intercept)", "x"))), w = structure(c(1, 1, 1, 0.5, 
-0.29999999999999999, 0.10000000000000001), dim = 3:2, dimnames = list(
    NULL, c("(Intercept)", "w"))), y = c(1, NA, 0), d = c(1, 
0, 1)), outer = 7L, inner = 13L, loglik = -3.3904979043503274, 
    score = c(0.75080292343635746, 0.71003447663193897, -0.080843405680819533, 
    -0.84171169465332374, -0.42274656946651085, -0.35813790793734523, 
    -0.015896632731514496, -0.072859239723542177)), `selection_probit 1` = list(
    theta = c(`selection.(Intercept)` = -0.14999999999999999, 
    selection.w = 0.29999999999999999, `outcome.(Intercept)` = 0.25, 
    outcome.x = -0.20000000000000001, sd_outcome = 0.40000000000000002, 
    sd_selection = 0.45000000000000001, cor_id = -0.34999999999999998, 
    cor_error = 0.25), panel = list(x = structure(c(1, 1, 1, 
    -0.40000000000000002, 0.29999999999999999, 0.69999999999999996
    ), dim = 3:2, dimnames = list(NULL, c("(Intercept)", "x"))), 
        w = structure(c(1, 1, 1, 0.5, -0.29999999999999999, 0.10000000000000001
        ), dim = 3:2, dimnames = list(NULL, c("(Intercept)", 
        "w"))), y = c(1, NA, 0), d = c(1, 0, 1)), outer = 7L, 
    inner = 13L, loglik = -3.4827571642055419, score = c(0.77098458560102212, 
    0.63463811864787867, -0.37382775570838178, -0.90963015072733944, 
    -0.31485104124205215, -0.29906625311083224, -0.046612948692259404, 
    -0.39851484987232338)), `selection_poisson -1` = list(theta = c(`selection.(Intercept)` = -0.14999999999999999, 
selection.w = 0.29999999999999999, `outcome.(Intercept)` = 0.25, 
outcome.x = -0.20000000000000001, sd_outcome = 0.40000000000000002, 
sd_selection = 0.45000000000000001, cor_id = 0.34999999999999998
), panel = list(x = structure(c(1, 1, 1, -0.40000000000000002, 
0.29999999999999999, 0.69999999999999996), dim = 3:2, dimnames = list(
    NULL, c("(Intercept)", "x"))), w = structure(c(1, 1, 1, 0.5, 
-0.29999999999999999, 0.10000000000000001), dim = 3:2, dimnames = list(
    NULL, c("(Intercept)", "w"))), y = c(1, NA, 2), d = c(1, 
0, 1)), outer = 7L, inner = 13L, loglik = -4.8884900466561092, 
    score = c(0.74061373308166389, 0.66594512059765965, 0.1578923165129735, 
    0.74494118097007911, -0.72691277109962793, -0.3194106025429182, 
    0.050525558237717372)), `selection_poisson 1` = list(theta = c(`selection.(Intercept)` = -0.14999999999999999, 
selection.w = 0.29999999999999999, `outcome.(Intercept)` = 0.25, 
outcome.x = -0.20000000000000001, sd_outcome = 0.40000000000000002, 
sd_selection = 0.45000000000000001, cor_id = -0.34999999999999998
), panel = list(x = structure(c(1, 1, 1, -0.40000000000000002, 
0.29999999999999999, 0.69999999999999996), dim = 3:2, dimnames = list(
    NULL, c("(Intercept)", "x"))), w = structure(c(1, 1, 1, 0.5, 
-0.29999999999999999, 0.10000000000000001), dim = 3:2, dimnames = list(
    NULL, c("(Intercept)", "w"))), y = c(1, NA, 2), d = c(1, 
0, 1)), outer = 7L, inner = 13L, loglik = -4.9121496386664498, 
    score = c(0.78237798304031758, 0.67000860828646636, 0.33910549137764123, 
    0.76120376617829899, -0.71683186628178375, -0.34619364587555973, 
    0.019497236284713306)), `selection_pln -1` = list(theta = c(`selection.(Intercept)` = -0.14999999999999999, 
selection.w = 0.29999999999999999, `outcome.(Intercept)` = 0.25, 
outcome.x = -0.20000000000000001, sd_outcome = 0.40000000000000002, 
sd_selection = 0.45000000000000001, cor_id = 0.34999999999999998, 
sd_error = 0.5, cor_error = -0.25), panel = list(x = structure(c(1, 
1, 1, -0.40000000000000002, 0.29999999999999999, 0.69999999999999996
), dim = 3:2, dimnames = list(NULL, c("(Intercept)", "x"))), 
    w = structure(c(1, 1, 1, 0.5, -0.29999999999999999, 0.10000000000000001
    ), dim = 3:2, dimnames = list(NULL, c("(Intercept)", "w"))), 
    y = c(1, NA, 2), d = c(1, 0, 1)), outer = 7L, inner = 13L, 
    loglik = -5.084783426786708, score = c(0.75012289674189314, 
    0.65478455898528298, 0.17645834023055454, 0.57007432167749017, 
    -0.57257366427543355, -0.32297629871367617, 0.031621613231063571, 
    -0.67399509954739434, 0.061954542220835122)), `selection_pln 1` = list(
    theta = c(`selection.(Intercept)` = -0.14999999999999999, 
    selection.w = 0.29999999999999999, `outcome.(Intercept)` = 0.25, 
    outcome.x = -0.20000000000000001, sd_outcome = 0.40000000000000002, 
    sd_selection = 0.45000000000000001, cor_id = -0.34999999999999998, 
    sd_error = 0.5, cor_error = 0.25), panel = list(x = structure(c(1, 
    1, 1, -0.40000000000000002, 0.29999999999999999, 0.69999999999999996
    ), dim = 3:2, dimnames = list(NULL, c("(Intercept)", "x"))), 
        w = structure(c(1, 1, 1, 0.5, -0.29999999999999999, 0.10000000000000001
        ), dim = 3:2, dimnames = list(NULL, c("(Intercept)", 
        "w"))), y = c(1, NA, 2), d = c(1, 0, 1)), outer = 7L, 
    inner = 13L, loglik = -5.066107987967813, score = c(0.75352754891791585, 
    0.68104909299161376, 0.030081564062103931, 0.53719016944819775, 
    -0.62617974746155425, -0.34173007531426092, -0.0040722984671557536, 
    -0.66373057904375521, 0.053960672406045143))), refined = list(
    probit_re = list(theta = c(`outcome.(Intercept)` = 0.25, 
    outcome.x = -0.20000000000000001, sd_outcome = 0.40000000000000002
    ), panel = list(x = structure(c(1, 1, 1, -0.40000000000000002, 
    0.29999999999999999, 0.69999999999999996), dim = 3:2, dimnames = list(
        NULL, c("(Intercept)", "x"))), w = structure(c(1, 1, 
    1, 0.5, -0.29999999999999999, 0.10000000000000001), dim = 3:2, dimnames = list(
        NULL, c("(Intercept)", "w"))), y = c(1, 0, 0), d = c(1, 
    1, 1)), outer = 31L, inner = 25L, loglik = -2.1473865321035404, 
        score = c(-0.92510669197616247, -1.0828343362830439, 
        -0.22831982091849284)), poisson_re = list(theta = c(`outcome.(Intercept)` = 0.25, 
    outcome.x = -0.20000000000000001, sd_outcome = 0.40000000000000002
    ), panel = list(x = structure(c(1, 1, 1, -0.40000000000000002, 
    0.29999999999999999, 0.69999999999999996), dim = 3:2, dimnames = list(
        NULL, c("(Intercept)", "x"))), w = structure(c(1, 1, 
    1, 0.5, -0.29999999999999999, 0.10000000000000001), dim = 3:2, dimnames = list(
        NULL, c("(Intercept)", "w"))), y = c(1, 0, 2), d = c(1, 
    1, 1)), outer = 31L, inner = 25L, loglik = -4.0555380134174159, 
        score = c(-0.56922028448762729, 0.43551007680608395, 
        -0.76868078952920471)), pln_re = list(theta = c(`outcome.(Intercept)` = 0.25, 
    outcome.x = -0.20000000000000001, sd_outcome = 0.40000000000000002, 
    sd_error = 0.5), panel = list(x = structure(c(1, 1, 1, -0.40000000000000002, 
    0.29999999999999999, 0.69999999999999996), dim = 3:2, dimnames = list(
        NULL, c("(Intercept)", "x"))), w = structure(c(1, 1, 
    1, 0.5, -0.29999999999999999, 0.10000000000000001), dim = 3:2, dimnames = list(
        NULL, c("(Intercept)", "w"))), y = c(1, 0, 2), d = c(1, 
    1, 1)), outer = 31L, inner = 25L, loglik = -4.1679846487288179, 
        score = c(-0.63650852243618317, 0.29177213281044395, 
        -0.59031122759562971, -0.47026491816174776)), selection_gaussian = list(
        theta = c(`selection.(Intercept)` = -0.14999999999999999, 
        selection.w = 0.29999999999999999, `outcome.(Intercept)` = 0.25, 
        outcome.x = -0.20000000000000001, sd_outcome = 0.40000000000000002, 
        sd_selection = 0.45000000000000001, cor_id = -0.34999999999999998, 
        sd_error = 0.5, cor_error = 0.25), panel = list(x = structure(c(1, 
        1, 1, -0.40000000000000002, 0.29999999999999999, 0.69999999999999996
        ), dim = 3:2, dimnames = list(NULL, c("(Intercept)", 
        "x"))), w = structure(c(1, 1, 1, 0.5, -0.29999999999999999, 
        0.10000000000000001), dim = 3:2, dimnames = list(NULL, 
            c("(Intercept)", "w"))), y = c(0.20000000000000001, 
        NA, -0.40000000000000002), d = c(1, 0, 1)), outer = 31L, 
        inner = 25L, loglik = -3.3235770294396403, score = c(0.71598811234992388, 
        0.65799445770828902, -1.2979422240656009, -1.0848782190816157, 
        -0.57817356106455675, -0.23070597944896326, -0.18586096485899176, 
        -2.0917399159658712, -0.4087276562497163)), selection_probit = list(
        theta = c(`selection.(Intercept)` = -0.14999999999999999, 
        selection.w = 0.29999999999999999, `outcome.(Intercept)` = 0.25, 
        outcome.x = -0.20000000000000001, sd_outcome = 0.40000000000000002, 
        sd_selection = 0.45000000000000001, cor_id = -0.34999999999999998, 
        cor_error = 0.25), panel = list(x = structure(c(1, 1, 
        1, -0.40000000000000002, 0.29999999999999999, 0.69999999999999996
        ), dim = 3:2, dimnames = list(NULL, c("(Intercept)", 
        "x"))), w = structure(c(1, 1, 1, 0.5, -0.29999999999999999, 
        0.10000000000000001), dim = 3:2, dimnames = list(NULL, 
            c("(Intercept)", "w"))), y = c(1, NA, 0), d = c(1, 
        0, 1)), outer = 31L, inner = 25L, loglik = -3.4827506390812766, 
        score = c(0.77094848285813067, 0.63463120421987762, -0.37381899819250058, 
        -0.90962776432271097, -0.31484134506761158, -0.29891965316027869, 
        -0.046630359363361656, -0.39849469509030866)), selection_poisson = list(
        theta = c(`selection.(Intercept)` = -0.14999999999999999, 
        selection.w = 0.29999999999999999, `outcome.(Intercept)` = 0.25, 
        outcome.x = -0.20000000000000001, sd_outcome = 0.40000000000000002, 
        sd_selection = 0.45000000000000001, cor_id = -0.34999999999999998
        ), panel = list(x = structure(c(1, 1, 1, -0.40000000000000002, 
        0.29999999999999999, 0.69999999999999996), dim = 3:2, dimnames = list(
            NULL, c("(Intercept)", "x"))), w = structure(c(1, 
        1, 1, 0.5, -0.29999999999999999, 0.10000000000000001), dim = 3:2, dimnames = list(
            NULL, c("(Intercept)", "w"))), y = c(1, NA, 2), d = c(1, 
        0, 1)), outer = 31L, inner = 25L, loglik = -4.9121445552480969, 
        score = c(0.78234674192805942, 0.67000473288694895, 0.33916907048303097, 
        0.76120947194786848, -0.71669006030799653, -0.3461566840616937, 
        0.01952804121459556)), selection_pln = list(theta = c(`selection.(Intercept)` = -0.14999999999999999, 
    selection.w = 0.29999999999999999, `outcome.(Intercept)` = 0.25, 
    outcome.x = -0.20000000000000001, sd_outcome = 0.40000000000000002, 
    sd_selection = 0.45000000000000001, cor_id = -0.34999999999999998, 
    sd_error = 0.5, cor_error = 0.25), panel = list(x = structure(c(1, 
    1, 1, -0.40000000000000002, 0.29999999999999999, 0.69999999999999996
    ), dim = 3:2, dimnames = list(NULL, c("(Intercept)", "x"))), 
        w = structure(c(1, 1, 1, 0.5, -0.29999999999999999, 0.10000000000000001
        ), dim = 3:2, dimnames = list(NULL, c("(Intercept)", 
        "w"))), y = c(1, NA, 2), d = c(1, 0, 1)), outer = 31L, 
        inner = 25L, loglik = -5.0661045568210792, score = c(0.75348505893484996, 
        0.68104320816366937, 0.030089360976077298, 0.53719136625467145, 
        -0.62616063261106736, -0.34166303638334422, -0.0040703308175334824, 
        -0.66373588334430522, 0.053961150561174703))))
