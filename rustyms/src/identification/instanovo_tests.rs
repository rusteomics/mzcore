#![allow(clippy::missing_panics_doc)]
use crate::identification::{test_format, InstaNovoData, InstaNovoVersion};
use std::io::BufReader;

#[test]
fn instanovo() {
    match test_format::<InstaNovoData>(
        BufReader::new(INSTANOVO_V1_0_0.as_bytes()),
        None,
        true,
        true,
        Some(InstaNovoVersion::V1_0_0),
    ) {
        Ok(n) => assert_eq!(n, 21),
        Err(e) => {
            println!("{e}");
            panic!("Failed identified peptides test");
        }
    }
}

const INSTANOVO_V1_0_0: &str = r#"scan_number,precursor_mz,precursor_charge,experiment_name,spectrum_id,preds,preds_tokenised,log_probs,token_log_probs
0,1353.116333007813,4,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp:0,LKVKVILEAEPS(+79.97)EEEEEEEEEEEEEEEEEEEEEEEEKEEK,"L, K, V, K, V, I, L, E, A, E, P, S(+79.97), E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, K, E, E, K",-47.14482498168945,"[-1.1898047924041748, -1.2532058954238892, -1.4706779718399048, -1.578391671180725, -1.910727858543396, -0.4288635551929474, -0.10262472927570343, -0.2335159033536911, -0.3816143870353699, -0.1399289071559906, -0.2679944634437561, -0.37487441301345825, -0.28003591299057007, -0.29957395792007446, -0.6062297224998474, -1.0798466205596924, -1.3055310249328613, -1.1969765424728394, -0.8466325402259827, -0.7559331655502319, -0.8520379066467285, -1.1635522842407227, -1.5230286121368408, -1.5223480463027954, -1.3874539136886597, -1.3555835485458374, -1.3308098316192627, -1.461938738822937, -1.292738437652588, -1.7667877674102783, -1.8383617401123047, -1.924727439880371, -1.5695301294326782, -1.4049240350723267, -1.2322568893432617, -1.1730256080627441, -0.09055394679307938, -3.6145036220550537, -1.8250231742858887, -3.1126556396484375]"
1,1353.116333007813,4,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp:1,SSSSTSGS(+79.97)DC(+57.02)DGVHVEPEEEDMES(+79.97)EDEDEDEDLVTSTTSK,"S, S, S, S, T, S, G, S(+79.97), D, C(+57.02), D, G, V, H, V, E, P, E, E, E, D, M, E, S(+79.97), E, D, E, D, E, D, E, D, L, V, T, S, T, T, S, K",-33.68059539794922,"[-0.07847003638744354, -0.5382011532783508, -0.38900521397590637, -0.3034709692001343, -0.14283983409404755, -0.004337664693593979, -0.005992896854877472, -0.443778932094574, -2.0242862701416016, -1.3575999736785889, -1.3120659589767456, -0.8160025477409363, -0.9171149730682373, -0.1492014229297638, -0.28191035985946655, -1.03749680519104, -0.7952876687049866, -0.11219097673892975, -0.6492378115653992, -0.0394880585372448, -0.35166993737220764, -0.031147046014666557, -0.1014118641614914, -0.8919384479522705, -0.5123100876808167, -1.5009464025497437, -0.7995803952217102, -0.8618601560592651, -0.74873948097229, -1.0185350179672241, -1.0649776458740234, -1.9072270393371582, -1.3915926218032837, -1.454406499862671, -0.8135812878608704, -0.7613984942436218, -1.0519371032714844, -2.232295274734497, -2.410634994506836, -2.3764281272888184]"
2,1216.224487304688,2,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp:2,LSPS(+79.97)PPLLPPPPPPVPLPPLPPK,"L, S, P, S(+79.97), P, P, L, L, P, P, P, P, P, P, V, P, L, P, P, L, P, P, K",-234.313232421875,"[-0.48280513286590576, -1.9982898235321045, -1.145535945892334, -1.2481436729431152, -1.242499589920044, -0.4877724349498749, -0.8774506449699402, -0.7018945217132568, -0.9126332998275757, -1.0259754657745361, -1.1820683479309082, -0.46859246492385864, -1.7760021686553955, -0.9172687530517578, -0.9521878361701965, -0.9288878440856934, -0.02773796021938324, -0.20492245256900787, -1.173280954360962, -0.6005244255065918, -0.05487779155373573, -0.6778227686882019, -2.4466471672058105]"
3,1216.224487304688,2,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp:3,SAPS(+79.97)PELLDLPPLLPRS(+79.97)PLPK,"S, A, P, S(+79.97), P, E, L, L, D, L, P, P, L, L, P, R, S(+79.97), P, L, P, K",-255.28448486328125,"[-0.05731990188360214, -0.6153630018234253, -1.2121751308441162, -1.229188323020935, -1.457108736038208, -1.4648027420043945, -1.4639601707458496, -1.4428002834320068, -0.28625816106796265, -0.0845002830028534, -1.1018378734588623, -1.210679531097412, -2.060678482055664, -2.2136943340301514, -0.7935038805007935, -1.8666895627975464, -1.7751319408416748, -0.8543295860290527, -1.4944884777069092, -1.0473661422729492, -1.3965165615081787]"
4,1349.11279296875,4,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp:4,GGGGGGGGGGGGGGGGGGGGGGGGAGEGAGADRS(+79.97)PQRPGR,"G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, G, A, G, E, G, A, G, A, D, R, S(+79.97), P, Q, R, P, G, R",-34.533355712890625,"[-0.6063277721405029, -0.6077976226806641, -0.08624760806560516, -0.5748220086097717, -1.5442652702331543, -1.4100384712219238, -1.241212010383606, -1.0145214796066284, -1.129785418510437, -1.2704263925552368, -1.169306755065918, -1.5771288871765137, -1.5485100746154785, -1.499711036682129, -1.2069355249404907, -1.5308293104171753, -1.247301697731018, -1.4038337469100952, -1.5108823776245117, -1.336652398109436, -1.2145107984542847, -0.9875333309173584, -0.7432393431663513, -0.5624361634254456, -0.49025505781173706, -0.519905686378479, -0.6296592354774475, -0.7260276675224304, -0.6793420910835266, -0.5041154623031616, -0.3764101266860962, -0.33699753880500793, -0.4097467362880707, -0.4173401892185211, -0.432345449924469, -0.4027464687824249, -0.3778119683265686, -0.3786942958831787, -0.40598759055137634, -0.42171356081962585]"
5,1349.11279296875,4,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp:5,STSTFQALQDELISLHPPEEEEEEEEEEES(+79.97)EDEEVTKEEK,"S, T, S, T, F, Q, A, L, Q, D, E, L, I, S, L, H, P, P, E, E, E, E, E, E, E, E, E, E, E, S(+79.97), E, D, E, E, V, T, K, E, E, K",-38.61011505126953,"[-0.27568867802619934, -1.4995495080947876, -1.4445313215255737, -1.3065897226333618, -0.9861798882484436, -0.8434268236160278, -0.7076863050460815, -1.287238359451294, -0.5336244106292725, -0.8561721444129944, -1.3078498840332031, -0.623564600944519, -0.7130038738250732, -0.8667600750923157, -0.8819077610969543, -1.0462291240692139, -1.1810555458068848, -0.9374122023582458, -0.817981481552124, -0.6692118644714355, -0.6896576285362244, -1.285140872001648, -1.6251142024993896, -0.3567267954349518, -0.46701139211654663, -1.557365894317627, -0.712128221988678, -1.579364538192749, -0.18570637702941895, -1.0566295385360718, -0.5056825876235962, -1.7716294527053833, -0.8920828700065613, -0.27432313561439514, -2.4893293380737305, -0.9152501225471497, -1.13835871219635, -0.7423431873321533, -1.0278247594833374, -0.552785336971283]"
6,1357.365478515625,4,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp:6,LSS(+79.97)EVTSSLSTSMGLHDEDDDFFDRDEDLKEDLVTISGTS,"L, S, S(+79.97), E, V, T, S, S, L, S, T, S, M, G, L, H, D, E, D, D, D, F, F, D, R, D, E, D, L, K, E, D, L, V, T, I, S, G, T, S",-55.809932708740234,"[-1.0700087547302246, -0.45235228538513184, -0.13142752647399902, -0.034827232360839844, -0.33834248781204224, -0.594786524772644, -0.6601390242576599, -0.8261093497276306, -1.635099172592163, -1.7555067539215088, -0.9311437606811523, -1.4158742427825928, -1.4413493871688843, -1.2220157384872437, -1.0625641345977783, -1.1492520570755005, -0.9419587850570679, -1.2509950399398804, -1.2027850151062012, -1.6362003087997437, -1.3815306425094604, -0.997757613658905, -0.7760587334632874, -1.0876846313476562, -1.5182546377182007, -1.0099989175796509, -0.26469892263412476, -1.6592745780944824, -1.8463445901870728, -0.8798919916152954, -1.8809072971343994, -0.30079385638237, -4.87946891784668, -0.3360159695148468, -0.27439162135124207, -1.700141429901123, -4.711390972137451, -5.989498615264893, -2.9504804611206055, -1.6126152276992798]"
7,1357.365478515625,4,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp:7,EEEEEEEEDEEEEEEEEEEEEEEEEEDERQS(+79.97)PVVLDSSSK,"E, E, E, E, E, E, E, E, D, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, D, E, R, Q, S(+79.97), P, V, V, L, D, S, S, S, K",-38.34630584716797,"[-0.15063245594501495, -0.8325968384742737, -0.9667963981628418, -0.26966795325279236, -1.1486272811889648, -1.0617644786834717, -0.9288771152496338, -0.10934460163116455, -0.510221540927887, -0.7267850041389465, -0.8313274383544922, -0.8427366018295288, -1.0229275226593018, -0.8890469074249268, -0.29960405826568604, -1.4826282262802124, -0.5701924562454224, -1.3117291927337646, -0.8567087650299072, -0.458834171295166, -0.30649030208587646, -0.5777329206466675, -1.1479213237762451, -1.302885890007019, -1.3883767127990723, -1.2489534616470337, -0.8490742444992065, -0.8153007626533508, -0.9253391623497009, -1.1632267236709595, -1.2933614253997803, -1.2381542921066284, -1.1385905742645264, -0.8864232897758484, -1.0033221244812012, -1.4086196422576904, -1.6194214820861816, -1.6751060485839844, -1.643905758857727, -1.4430569410324097]"
8,1349.367553710938,4,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp:8,SASISWS(+79.97)PVKC(+57.02)SVWDPDHDPDHDDADYWMKDDLGLTSGTS,"S, A, S, I, S, W, S(+79.97), P, V, K, C(+57.02), S, V, W, D, P, D, H, D, P, D, H, D, D, A, D, Y, W, M, K, D, D, L, G, L, T, S, G, T, S",-65.87045288085938,"[-0.406949520111084, -0.13111849129199982, -0.09339886158704758, -0.006496618967503309, -0.29331085085868835, -1.0819785594940186, -1.07078218460083, -1.3674098253250122, -1.018458604812622, -1.596422553062439, -1.7888071537017822, -1.711216926574707, -1.337314486503601, -1.3132869005203247, -1.4704970121383667, -1.6973981857299805, -1.139493703842163, -0.4161772131919861, -1.207051396369934, -1.0346803665161133, -1.3333330154418945, -1.3219164609909058, -2.0319883823394775, -1.1902544498443604, -0.9279484152793884, -1.4660371541976929, -1.876696228981018, -0.5972193479537964, -0.32778477668762207, -0.5418118238449097, -7.632816314697266, -5.163928508758545, -1.9493813514709473, -0.7559415698051453, -2.515380382537842, -4.5958428382873535, -2.766427993774414, -3.0831775665283203, -2.661452293395996, -2.9488584995269775]"
9,1349.367553710938,4,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp:9,SSSSTSTSSSSTSSTDSPSGSEDRS(+79.97)PS(+79.97)PEKENEVTSETSK,"S, S, S, S, T, S, T, S, S, S, S, T, S, S, T, D, S, P, S, G, S, E, D, R, S(+79.97), P, S(+79.97), P, E, K, E, N, E, V, T, S, E, T, S, K",-47.47643280029297,"[-0.16404178738594055, -0.6305190920829773, -0.9443010091781616, -0.8775668144226074, -0.3844778537750244, -0.015940023586153984, -0.01927993819117546, -1.3933794498443604, -0.9012462496757507, -0.7598123550415039, -0.6868025660514832, -0.45880424976348877, -1.1288713216781616, -0.6652687788009644, -1.4917452335357666, -1.0591261386871338, -0.2025306075811386, -0.1081993505358696, -0.6764559149742126, -1.0689442157745361, -0.6504117250442505, -0.8548192977905273, -0.9655941128730774, -1.2193018198013306, -1.7102758884429932, -1.5224287509918213, -1.4552305936813354, -1.2006592750549316, -0.8853554725646973, -1.516101598739624, -1.3193846940994263, -1.0446220636367798, -1.2151552438735962, -0.8904541730880737, -3.5030975341796875, -2.7242209911346436, -3.4343440532684326, -3.248047113418579, -2.2567732334136963, -2.222841262817383]"
10,1353.370849609375,4,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp:10,VTSLSSSPLKPVLEAWIHDEEMDS(+79.97)EDEEEKWDDLYASSTS,"V, T, S, L, S, S, S, P, L, K, P, V, L, E, A, W, I, H, D, E, E, M, D, S(+79.97), E, D, E, E, E, K, W, D, D, L, Y, A, S, S, T, S",-44.6576042175293,"[-0.17519927024841309, -0.07632763683795929, -0.34875380992889404, -0.042356234043836594, -0.013340969569981098, -1.274790644645691, -0.7739483714103699, -1.4739679098129272, -1.3003861904144287, -1.310808777809143, -1.522390365600586, -1.3426114320755005, -1.4101788997650146, -0.8371820449829102, -0.7036910057067871, -0.3988977372646332, -0.9187065958976746, -0.23145386576652527, -1.8500487804412842, -0.5929505228996277, -0.8790079951286316, -1.0516234636306763, -0.8873531818389893, -1.32279634475708, -1.5349055528640747, -0.9580303430557251, -1.2929342985153198, -1.8876397609710693, -1.7891830205917358, -0.6723992824554443, -0.13665027916431427, -2.408869981765747, -0.2714077830314636, -0.340230792760849, -0.9830684661865234, -0.8312598466873169, -0.5362821221351624, -6.1171956062316895, -1.9116570949554443, -2.2471208572387695]"
11,1353.370849609375,4,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp:11,EEEEEEEEEEEEEEEEEEEEEEEEEEMERQS(+79.97)PVVLDSSSK,"E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, E, M, E, R, Q, S(+79.97), P, V, V, L, D, S, S, S, K",-42.54265213012695,"[-0.2640531063079834, -1.1067227125167847, -0.7537723779678345, -0.21282437443733215, -0.6444227695465088, -0.5025171637535095, -0.11948210000991821, -0.14578743278980255, -0.8067586421966553, -0.39255577325820923, -0.37441569566726685, -0.47786927223205566, -1.2774714231491089, -0.943406343460083, -1.488488793373108, -1.0032429695129395, -1.2166993618011475, -1.558054804801941, -1.3079780340194702, -1.1804583072662354, -0.8643602132797241, -0.947065532207489, -1.3344824314117432, -1.4143755435943604, -1.4135427474975586, -1.4270200729370117, -1.3108503818511963, -1.2149748802185059, -1.1543819904327393, -1.2545700073242188, -1.2425501346588135, -1.3000712394714355, -1.426724910736084, -1.437530279159546, -1.484959602355957, -1.4973523616790771, -1.476752758026123, -1.4833042621612549, -1.5526419878005981, -1.5281561613082886]"
12,1348.869384765625,4,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp:12,SSGGYWGSPGKYSSHHFMGYGYDDMYWKPSGSILVTSGTS,"S, S, G, G, Y, W, G, S, P, G, K, Y, S, S, H, H, F, M, G, Y, G, Y, D, D, M, Y, W, K, P, S, G, S, I, L, V, T, S, G, T, S",-71.5739974975586,"[-0.609157383441925, -0.18590812385082245, -0.18178929388523102, -0.07940386980772018, -0.2525795102119446, -0.5205252766609192, -0.7315086126327515, -1.489088535308838, -1.013169527053833, -1.2799092531204224, -1.231743335723877, -1.1684287786483765, -0.789280116558075, -1.224703311920166, -0.8068684339523315, -2.057281494140625, -1.8166402578353882, -1.4321420192718506, -1.009615182876587, -1.355675458908081, -1.7395979166030884, -0.051246218383312225, -1.2518012523651123, -1.2463343143463135, -1.091395616531372, -0.7693159580230713, -1.219667911529541, -0.8084853887557983, -0.7904479503631592, -1.9142206907272339, -5.524913787841797, -4.524893760681152, -3.785367727279663, -2.1374785900115967, -2.288924217224121, -1.1607757806777954, -7.280113697052002, -5.7617974281311035, -5.605969429016113, -3.385830879211426]"
13,1348.869384765625,4,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp:13,SSSSSNS(+79.97)QPPPLLEGERS(+79.97)EEEENGDGHSLFGNEVTSTSEK,"S, S, S, S, S, N, S(+79.97), Q, P, P, P, L, L, E, G, E, R, S(+79.97), E, E, E, E, N, G, D, G, H, S, L, F, G, N, E, V, T, S, T, S, E, K",-54.074180603027344,"[-0.20559526979923248, -1.5363306999206543, -0.9894521236419678, -0.7482156753540039, -0.15950946509838104, -0.0037842821329832077, -0.00392345804721117, -1.0287175178527832, -0.5543951392173767, -1.23534095287323, -0.43195152282714844, -1.4651038646697998, -2.116406202316284, -1.0090972185134888, -1.236779808998108, -1.435746431350708, -1.3745900392532349, -1.1407020092010498, -1.6297179460525513, -1.5244817733764648, -0.9910011291503906, -1.4753377437591553, -1.227499008178711, -1.6027615070343018, -1.2512940168380737, -1.066943883895874, -1.301374912261963, -1.7629523277282715, -1.4931827783584595, -1.3698959350585938, -1.8839346170425415, -1.6111435890197754, -1.7559510469436646, -1.1649703979492188, -1.7128162384033203, -1.4801974296569824, -2.850574016571045, -3.051332473754883, -2.797579765319824, -2.3935980796813965]"
14,1349.122436523438,4,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp:14,LSDFSISFWHS(+79.97)EEEDDDEEEDGGKKPHSYSFDSLYASSTS,"L, S, D, F, S, I, S, F, W, H, S(+79.97), E, E, E, D, D, D, E, E, E, D, G, G, K, K, P, H, S, Y, S, F, D, S, L, Y, A, S, S, T, S",-48.55881118774414,"[-0.37943023443222046, -0.05031585320830345, -0.4418899416923523, -0.5165826678276062, -0.329465389251709, -0.634541928768158, -0.3490736186504364, -1.7085533142089844, -1.7917423248291016, -0.2783038318157196, -0.9645243287086487, -0.9889963865280151, -0.803525984287262, -1.2197767496109009, -0.748043954372406, -0.9806399345397949, -2.236006259918213, -0.8217480182647705, -1.6584222316741943, -0.941186249256134, -0.5847094058990479, -1.0577900409698486, -0.9975808262825012, -0.7513730525970459, -1.3372840881347656, -1.5886975526809692, -0.7404971122741699, -1.5333542823791504, -0.4998478293418884, -1.0035555362701416, -1.750340461730957, -1.3353850841522217, -1.6656231880187988, -0.9581365585327148, -1.0396952629089355, -5.255682468414307, -1.6375943422317505, -3.4260737895965576, -2.195528030395508, -1.3572865724563599]"
15,1349.122436523438,4,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp:15,TSTSSTSTEAQVHAEAQAHPEDRS(+79.97)PEEEAEGAEEEKEAEK,"T, S, T, S, S, T, S, T, E, A, Q, V, H, A, E, A, Q, A, H, P, E, D, R, S(+79.97), P, E, E, E, A, E, G, A, E, E, E, K, E, A, E, K",-35.47725296020508,"[-0.12350436300039291, -0.22369863092899323, -0.04920090734958649, -1.4021598100662231, -0.7673075795173645, -0.6899320483207703, -0.23632514476776123, -0.21554967761039734, -0.5195893049240112, -0.2675805985927582, -0.20016951858997345, -0.025143157690763474, -0.34623104333877563, -1.402597427368164, -0.7401972413063049, -0.7556778788566589, -1.0256924629211426, -0.9343301057815552, -1.0677529573440552, -0.9464545249938965, -1.5589196681976318, -1.0827124118804932, -0.7851864695549011, -1.941973328590393, -1.2625460624694824, -0.7867250442504883, -1.538809061050415, -1.6573245525360107, -1.3993407487869263, -1.3885537385940552, -1.1303845643997192, -1.4581077098846436, -1.1827868223190308, -0.11755254119634628, -0.853971004486084, -0.8964711427688599, -1.5724472999572754, -0.8693927526473999, -1.1797130107879639, -0.8752410411834717]"
16,1353.110473632813,4,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp:16,LSSTSGS(+79.97)MSNWPPPPDEEGEKGEKPEDDAWLSLSSASGTS,"L, S, S, T, S, G, S(+79.97), M, S, N, W, P, P, P, P, D, E, E, G, E, K, G, E, K, P, E, D, D, A, W, L, S, L, S, S, A, S, G, T, S",-65.9991683959961,"[-0.4931744933128357, -0.10294269025325775, -0.6120452880859375, -0.0013291343348100781, -0.15360486507415771, -0.5557436943054199, -0.325517475605011, -0.10480454564094543, -0.7972138524055481, -0.7757546305656433, -1.6229450702667236, -1.4787508249282837, -1.683245062828064, -1.440244197845459, -1.890653133392334, -0.7446271181106567, -0.7425130009651184, -1.9224246740341187, -1.027151346206665, -1.1922221183776855, -1.9122194051742554, -0.3947625756263733, -1.3419197797775269, -1.0347391366958618, -1.286238431930542, -1.3018429279327393, -0.17000213265419006, -0.5778448581695557, -1.8594253063201904, -2.606919288635254, -2.1031293869018555, -2.234875440597534, -1.6651183366775513, -7.347693920135498, -3.3774211406707764, -2.649264335632324, -3.1567208766937256, -4.120492458343506, -3.9556376934051514, -5.235989570617676]"
17,1353.110473632813,4,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp:17,DQAAGQDDLNPPLRTGDLGIPPNPEDRS(+79.97)PS(+79.97)PEPIYNSEGK,"D, Q, A, A, G, Q, D, D, L, N, P, P, L, R, T, G, D, L, G, I, P, P, N, P, E, D, R, S(+79.97), P, S(+79.97), P, E, P, I, Y, N, S, E, G, K",-34.3138313293457,"[-1.014553427696228, -1.1708364486694336, -0.13190269470214844, -0.14145584404468536, -0.19182397425174713, -0.007062352728098631, -0.00012492353562265635, -0.0034641751553863287, -0.0003505330823827535, -4.7444173105759546e-05, -0.00016926287207752466, -1.311301275563892e-06, -0.0008216104470193386, -0.0008131535141728818, -1.6689286894688848e-06, -1.4305104514278355e-06, -8.344646857949556e-07, -9.77468371274881e-05, -8.928377064876258e-05, -3.9457496313843876e-05, -0.001327943871729076, -7.128461584215984e-05, -3.731181277544238e-05, -2.8132995794294402e-05, -7.748573807475623e-06, -0.010917454957962036, -0.7962529063224792, -0.012375117279589176, -1.0862914323806763, -0.34013527631759644, -0.16126516461372375, -0.42035648226737976, -0.5432164669036865, -1.0152329206466675, -0.0848957896232605, -9.32642650604248, -3.798736810684204, -5.69951868057251, -4.028213977813721, -4.324866771697998]"
18,1353.363037109375,4,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp:18,GSGSEGGKPGWGWTHSRGYSFPDYDSNHSSSYVTVTSGST,"G, S, G, S, E, G, G, K, P, G, W, G, W, T, H, S, R, G, Y, S, F, P, D, Y, D, S, N, H, S, S, S, Y, V, T, V, T, S, G, S, T",-57.210357666015625,"[-0.08833134174346924, -0.010024912655353546, -0.000341476290486753, -0.0877658873796463, -0.018817594274878502, -1.3214493989944458, -0.3156588673591614, -0.9006086587905884, -1.6607422828674316, -1.2461398839950562, -1.5090794563293457, -0.9961817264556885, -1.3333173990249634, -1.5122567415237427, -0.4001505672931671, -1.5853387117385864, -1.7091730833053589, -0.7536590099334717, -1.4411301612854004, -1.387695550918579, -0.7693523168563843, -0.2747490704059601, -0.81729656457901, -0.8858135342597961, -1.272748589515686, -0.744476318359375, -2.2402219772338867, -1.3256865739822388, -1.1227340698242188, -1.442939043045044, -3.3821933269500732, -2.627223014831543, -2.360701322555542, -0.5053297281265259, -0.08112901449203491, -4.031546115875244, -3.548017978668213, -3.7473630905151367, -4.817002296447754, -2.935967206954956]"
19,1353.363037109375,4,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp:19,GGGGSAGGGGGAAGGRPS(+79.97)PPQENGKEPC(+57.02)EPSQS(+79.97)VTSTSEK,"G, G, G, G, S, A, G, G, G, G, G, A, A, G, G, R, P, S(+79.97), P, P, Q, E, N, G, K, E, P, C(+57.02), E, P, S, Q, S(+79.97), V, T, S, T, S, E, K",-34.095943450927734,"[-0.3195824921131134, -1.3936655521392822, -1.1236826181411743, -0.7992289066314697, -0.2012583613395691, -0.005753267090767622, -0.07621022313833237, -1.300059199333191, -1.4463657140731812, -1.1221880912780762, -0.6981537938117981, -1.6815447807312012, -0.3504905104637146, -1.444809913635254, -1.5379170179367065, -0.9235746264457703, -0.9498796463012695, -0.26199230551719666, -0.04976364225149155, -0.6369405388832092, -0.9404551982879639, -0.29938265681266785, -1.3071815967559814, -0.10240114480257034, -0.7078803777694702, -0.22844867408275604, -0.8289092779159546, -0.9006180763244629, -1.3896540403366089, -0.18864311277866364, -0.2939540445804596, -0.9182708859443665, -0.897458553314209, -1.0683497190475464, -0.998686671257019, -1.3441799879074097, -1.3926314115524292, -1.3044145107269287, -1.4491615295410156, -1.2121996879577637]"
20,1348.856201171875,4,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp,20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp:20,GSGSGRGS(+79.97)WGSGGGHSSYSQHFNKPIFLKPGQWTWISGTS,"G, S, G, S, G, R, G, S(+79.97), W, G, S, G, G, G, H, S, S, Y, S, Q, H, F, N, K, P, I, F, L, K, P, G, Q, W, T, W, I, S, G, T, S",-62.733360290527344,"[-0.3141501843929291, -0.14591962099075317, -0.4127911925315857, -0.0423421785235405, -0.1893325001001358, -1.3861042261123657, -0.8396965861320496, -1.3329275846481323, -1.982820987701416, -1.4959850311279297, -1.1631689071655273, -0.4160800278186798, -1.9834104776382446, -1.8969424962997437, -1.2355120182037354, -0.7272624373435974, -0.49926871061325073, -1.3329787254333496, -1.0549516677856445, -0.562333345413208, -0.6821218729019165, -0.908167839050293, -0.9920955896377563, -1.0761539936065674, -0.5029982328414917, -1.943763017654419, -1.5800557136535645, -1.903761863708496, -0.6185683608055115, -1.488539218902588, -0.8119556307792664, -3.762935161590576, -4.434981822967529, -2.770512104034424, -1.9104329347610474, -3.1854987144470215, -3.1550867557525635, -4.076790809631348, -4.346488952636719, -3.5684683322906494]""#;
