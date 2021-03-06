#include <gsl/gsl_sort_vector.h>
#include <math.h>
#include <munit.h>

#include <qdm.h>

#include "test.h"

static
MunitResult
test_qdm_ispline_matrix(
    MUNIT_UNUSED const MunitParameter params[],
    MUNIT_UNUSED void* fixture
)
{
  double taus_data[] = {
      0.010000000000000000,
      0.014999999999999999,
      0.020000000000000000,
      0.025000000000000001,
      0.029999999999999999,
      0.035000000000000003,
      0.040000000000000001,
      0.045000000000000005,
      0.050000000000000003,
      0.055000000000000000,
      0.060000000000000005,
      0.065000000000000002,
      0.069999999999999993,
      0.074999999999999997,
      0.080000000000000002,
      0.084999999999999992,
      0.089999999999999997,
      0.095000000000000001,
      0.099999999999999992,
      0.104999999999999996,
      0.110000000000000001,
      0.114999999999999991,
      0.119999999999999996,
      0.125000000000000000,
      0.130000000000000004,
      0.135000000000000009,
      0.140000000000000013,
      0.145000000000000018,
      0.150000000000000022,
      0.154999999999999999,
      0.160000000000000003,
      0.165000000000000008,
      0.170000000000000012,
      0.175000000000000017,
      0.180000000000000021,
      0.185000000000000026,
      0.190000000000000002,
      0.195000000000000007,
      0.200000000000000011,
      0.205000000000000016,
      0.210000000000000020,
      0.215000000000000024,
      0.220000000000000001,
      0.225000000000000006,
      0.230000000000000010,
      0.235000000000000014,
      0.240000000000000019,
      0.245000000000000023,
      0.250000000000000000,
      0.255000000000000004,
      0.260000000000000009,
      0.265000000000000013,
      0.270000000000000018,
      0.275000000000000022,
      0.280000000000000027,
      0.285000000000000031,
      0.290000000000000036,
      0.295000000000000040,
      0.299999999999999989,
      0.304999999999999993,
      0.309999999999999998,
      0.315000000000000002,
      0.320000000000000007,
      0.325000000000000011,
      0.330000000000000016,
      0.335000000000000020,
      0.340000000000000024,
      0.345000000000000029,
      0.350000000000000033,
      0.355000000000000038,
      0.360000000000000042,
      0.364999999999999991,
      0.369999999999999996,
      0.375000000000000000,
      0.380000000000000004,
      0.385000000000000009,
      0.390000000000000013,
      0.395000000000000018,
      0.400000000000000022,
      0.405000000000000027,
      0.410000000000000031,
      0.415000000000000036,
      0.420000000000000040,
      0.425000000000000044,
      0.429999999999999993,
      0.434999999999999998,
      0.440000000000000002,
      0.445000000000000007,
      0.450000000000000011,
      0.455000000000000016,
      0.460000000000000020,
      0.465000000000000024,
      0.470000000000000029,
      0.475000000000000033,
      0.480000000000000038,
      0.485000000000000042,
      0.489999999999999991,
      0.494999999999999996,
      0.500000000000000000,
      0.505000000000000004,
      0.510000000000000009,
      0.515000000000000013,
      0.520000000000000018,
      0.525000000000000022,
      0.530000000000000027,
      0.535000000000000031,
      0.540000000000000036,
      0.545000000000000040,
      0.550000000000000044,
      0.555000000000000049,
      0.560000000000000053,
      0.565000000000000058,
      0.570000000000000062,
      0.575000000000000067,
      0.580000000000000071,
      0.585000000000000075,
      0.589999999999999969,
      0.594999999999999973,
      0.599999999999999978,
      0.604999999999999982,
      0.609999999999999987,
      0.614999999999999991,
      0.619999999999999996,
      0.625000000000000000,
      0.630000000000000004,
      0.635000000000000009,
      0.640000000000000013,
      0.645000000000000018,
      0.650000000000000022,
      0.655000000000000027,
      0.660000000000000031,
      0.665000000000000036,
      0.670000000000000040,
      0.675000000000000044,
      0.680000000000000049,
      0.685000000000000053,
      0.690000000000000058,
      0.695000000000000062,
      0.700000000000000067,
      0.705000000000000071,
      0.710000000000000075,
      0.714999999999999969,
      0.719999999999999973,
      0.724999999999999978,
      0.729999999999999982,
      0.734999999999999987,
      0.739999999999999991,
      0.744999999999999996,
      0.750000000000000000,
      0.755000000000000004,
      0.760000000000000009,
      0.765000000000000013,
      0.770000000000000018,
      0.775000000000000022,
      0.780000000000000027,
      0.785000000000000031,
      0.790000000000000036,
      0.795000000000000040,
      0.800000000000000044,
      0.805000000000000049,
      0.810000000000000053,
      0.815000000000000058,
      0.820000000000000062,
      0.825000000000000067,
      0.830000000000000071,
      0.835000000000000075,
      0.840000000000000080,
      0.844999999999999973,
      0.849999999999999978,
      0.854999999999999982,
      0.859999999999999987,
      0.864999999999999991,
      0.869999999999999996,
      0.875000000000000000,
      0.880000000000000004,
      0.885000000000000009,
      0.890000000000000013,
      0.895000000000000018,
      0.900000000000000022,
      0.905000000000000027,
      0.910000000000000031,
      0.915000000000000036,
      0.920000000000000040,
      0.925000000000000044,
      0.930000000000000049,
      0.935000000000000053,
      0.940000000000000058,
      0.945000000000000062,
      0.950000000000000067,
      0.955000000000000071,
      0.960000000000000075,
      0.965000000000000080,
      0.969999999999999973,
      0.974999999999999978,
      0.979999999999999982,
      0.984999999999999987,
      0.989999999999999991,
  };

  gsl_vector_view taus = gsl_vector_view_array(taus_data, sizeof(taus_data) / sizeof(double));

  size_t spline_df = 3;

  double knots_data[] = {
      0, 0, 0,
      0.31526315789473686,
      0.47947368421052627,
      1, 1, 1,
  };

  gsl_vector_view knots = gsl_vector_view_array(knots_data, sizeof(knots_data) / sizeof(double));

  double ix_data[] = {
      1 , 0.0921721253598851 , 0.00194986550917209 , 6.61548977531158e-06 , 0                    , 0                    ,
      1 , 0.136054236987845  , 0.00434806829428824 , 2.23272779916766e-05 , 0                    , 0                    ,
      1 , 0.178498961847979  , 0.007660336343003   , 5.29239182024926e-05 , 0                    , 0                    ,
      1 , 0.219530235390408  , 0.0118605835877506  , 0.000103367027739243 , 0                    , 0                    ,
      1 , 0.259171993065254  , 0.0169227239609648  , 0.000178618223933413 , 0                    , 0                    ,
      1 , 0.297448170322639  , 0.0228206713950798  , 0.000283639124116484 , 0                    , 0                    ,
      1 , 0.334382702612684  , 0.0295283398225297  , 0.000423391345619941 , 0                    , 0                    ,
      1 , 0.36999952538551   , 0.0370196431757479  , 0.000602836505775268 , 0                    , 0                    ,
      1 , 0.40432257409124   , 0.0452684953871691  , 0.000826936221913948 , 0                    , 0                    ,
      1 , 0.437375784179994  , 0.0542488103892265  , 0.00110065211136746  , 0                    , 0                    ,
      1 , 0.469183091101894  , 0.0639345021143549  , 0.0014289457914673   , 0                    , 0                    ,
      1 , 0.499768430307062  , 0.074299484494988   , 0.00181677887954494  , 0                    , 0                    ,
      1 , 0.52915573724562   , 0.0853176714635595  , 0.00226911299293187  , 0                    , 0                    ,
      1 , 0.557368947367689  , 0.096962976952504   , 0.00279090974895957  , 0                    , 0                    ,
      1 , 0.58443199612339   , 0.109209314894255   , 0.00338713076495953  , 0                    , 0                    ,
      1 , 0.610368818962845  , 0.122030599221247   , 0.00406273765826322  , 0                    , 0                    ,
      1 , 0.635203351336176  , 0.135400743865913   , 0.00482269204620214  , 0                    , 0                    ,
      1 , 0.658959528693504  , 0.149293662760688   , 0.00567195554610777  , 0                    , 0                    ,
      1 , 0.681661286484951  , 0.163683269838005   , 0.00661548977531158  , 0                    , 0                    ,
      1 , 0.703332560160638  , 0.178543479030299   , 0.00765825635114506  , 0                    , 0                    ,
      1 , 0.723997285170687  , 0.193848204270003   , 0.00880521689093971  , 0                    , 0                    ,
      1 , 0.743679396965219  , 0.209571359489553   , 0.010061333012027    , 0                    , 0                    ,
      1 , 0.762402830994356  , 0.22568685862138    , 0.0114315663317384   , 0                    , 0                    ,
      1 , 0.78019152270822   , 0.242168615597921   , 0.0129208784674054   , 0                    , 0                    ,
      1 , 0.797069407556932  , 0.258990544351608   , 0.0145342310363595   , 0                    , 0                    ,
      1 , 0.813060420990613  , 0.276126558814875   , 0.0162765856559322   , 0                    , 0                    ,
      1 , 0.828188498459386  , 0.293550572920157   , 0.018152903943455    , 0                    , 0                    ,
      1 , 0.842477575413372  , 0.311236500599888   , 0.0201681475162593   , 0                    , 0                    ,
      1 , 0.855951587302691  , 0.329158255786501   , 0.0223272779916766   , 0                    , 0                    ,
      1 , 0.868634469577467  , 0.34728975241243    , 0.0246352569870384   , 0                    , 0                    ,
      1 , 0.88055015768782   , 0.36560490441011    , 0.0270970461196762   , 0                    , 0                    ,
      1 , 0.891722587083872  , 0.384077625711975   , 0.0297176070069215   , 0                    , 0                    ,
      1 , 0.902175693215745  , 0.402681830250458   , 0.0325019012661058   , 0                    , 0                    ,
      1 , 0.91193341153356   , 0.421391431957994   , 0.0354548905145605   , 0                    , 0                    ,
      1 , 0.921019677487438  , 0.440180344767016   , 0.0385815363696171   , 0                    , 0                    ,
      1 , 0.929458426527502  , 0.459022482609959   , 0.0418868004486072   , 0                    , 0                    ,
      1 , 0.937273594103872  , 0.477891759419257   , 0.0453756443688621   , 0                    , 0                    ,
      1 , 0.944489115666671  , 0.496762089127342   , 0.0490530297477134   , 0                    , 0                    ,
      1 , 0.95112892666602   , 0.51560738566665    , 0.0529239182024926   , 0                    , 0                    ,
      1 , 0.95721696255204   , 0.534401562969615   , 0.0569932713505312   , 0                    , 0                    ,
      1 , 0.962777158774853  , 0.55311853496867    , 0.0612660508091605   , 0                    , 0                    ,
      1 , 0.96783345078458   , 0.571732215596249   , 0.0657472181957122   , 0                    , 0                    ,
      1 , 0.972409774031344  , 0.590216518784787   , 0.0704417351275177   , 0                    , 0                    ,
      1 , 0.976530063965266  , 0.608545358466718   , 0.0753545632219085   , 0                    , 0                    ,
      1 , 0.980218256036467  , 0.626692648574474   , 0.080490664096216    , 0                    , 0                    ,
      1 , 0.983498285695068  , 0.644632303040492   , 0.0858549993677718   , 0                    , 0                    ,
      1 , 0.986394088391192  , 0.662338235797203   , 0.0914525306539073   , 0                    , 0                    ,
      1 , 0.98892959957496   , 0.679784360777043   , 0.097288219571954    , 0                    , 0                    ,
      1 , 0.991128754696493  , 0.696944591912446   , 0.103367027739243    , 0                    , 0                    ,
      1 , 0.993015489205913  , 0.713792843135845   , 0.109693916773107    , 0                    , 0                    ,
      1 , 0.994613738553342  , 0.730303028379674   , 0.116273848290876    , 0                    , 0                    ,
      1 , 0.995947438188901  , 0.746449061576367   , 0.123111783909883    , 0                    , 0                    ,
      1 , 0.997040523562712  , 0.76220485665836    , 0.130212685247458    , 0                    , 0                    ,
      1 , 0.997916930124896  , 0.777544327558084   , 0.137581513920933    , 0                    , 0                    ,
      1 , 0.998600593325575  , 0.792441388207974   , 0.14522323154764     , 0                    , 0                    ,
      1 , 0.99911544861487   , 0.806869952540466   , 0.15314279974491     , 0                    , 0                    ,
      1 , 0.999485431442904  , 0.820803934487991   , 0.161345180130074    , 0                    , 0                    ,
      1 , 0.999734477259796  , 0.834217247982985   , 0.169835334320465    , 0                    , 0                    ,
      1 , 0.99988652151567   , 0.847083806957881   , 0.178618223933413    , 0                    , 0                    ,
      1 , 0.999965499660646  , 0.859377525345114   , 0.18769881058625     , 0                    , 0                    ,
      1 , 0.999995347144847  , 0.871072317077117   , 0.197082055896307    , 0                    , 0                    ,
      1 , 0.999999999418393  , 0.882142096086324   , 0.206772921480917    , 0                    , 0                    ,
      1 , 1                  , 0.892567288378607   , 0.216773370700693    , 1.3804425623868e-06  , 0                    ,
      1 , 1                  , 0.902358831346444   , 0.227071319014193    , 1.19896359932405e-05 , 0                    ,
      1 , 1                  , 0.911536590801244   , 0.237650571101383    , 4.15685529897326e-05 , 0                    ,
      1 , 1                  , 0.920120433671027   , 0.248494931128135    , 9.98584029504342e-05 , 0                    ,
      1 , 1                  , 0.928130226883813   , 0.259588203260312    , 0.000196600395273916 , 0                    ,
      1 , 1                  , 0.935585837367622   , 0.270914191663784    , 0.00034153573935875  , 0                    ,
      1 , 1                  , 0.942507132050474   , 0.282456700504415    , 0.000544405644603507 , 0                    ,
      1 , 1                  , 0.948913977860389   , 0.294199533948074    , 0.000814951320406757 , 0                    ,
      1 , 1                  , 0.954826241725386   , 0.30612649616062     , 0.00116291397616707  , 0                    ,
      1 , 1                  , 0.960263790573486   , 0.318221391307929    , 0.00159803482128302  , 0                    ,
      1 , 1                  , 0.965246491332709   , 0.330468023555862    , 0.00213005506515318  , 0                    ,
      1 , 1                  , 0.969794210931074   , 0.342850197070287    , 0.00276871591717611  , 0                    ,
      1 , 1                  , 0.973926816296602   , 0.355351716017071    , 0.0035237585867504   , 0                    ,
      1 , 1                  , 0.977664174357312   , 0.367956384562077    , 0.0044049242832746   , 0                    ,
      1 , 1                  , 0.981026152041225   , 0.380648006871176    , 0.0054219542161473   , 0                    ,
      1 , 1                  , 0.984032616276359   , 0.393410387110234    , 0.00658458959476706  , 0                    ,
      1 , 1                  , 0.986703433990736   , 0.406227329445115    , 0.00790257162853245  , 0                    ,
      1 , 1                  , 0.989058472112375   , 0.419082638041686    , 0.00938564152684205  , 0                    ,
      1 , 1                  , 0.991117597569296   , 0.431960117065814    , 0.0110435404990944   , 0                    ,
      1 , 1                  , 0.992900677289519   , 0.444843570683369    , 0.0128860097546881   , 0                    ,
      1 , 1                  , 0.994427578201064   , 0.457716803060211    , 0.0149227905030218   , 0                    ,
      1 , 1                  , 0.995718167231951   , 0.470563618362209    , 0.0171636239534939   , 0                    ,
      1 , 1                  , 0.996792311310199   , 0.483367820755232    , 0.0196182513155031   , 0                    ,
      1 , 1                  , 0.997669877363829   , 0.496113214405141    , 0.0222964137984479   , 0                    ,
      1 , 1                  , 0.998370732320861   , 0.50878360347781     , 0.0252078526117269   , 0                    ,
      1 , 1                  , 0.998914743109315   , 0.5213627921391      , 0.0283623089647387   , 0                    ,
      1 , 1                  , 0.99932177665721    , 0.533834584554879    , 0.0317695240668818   , 0                    ,
      1 , 1                  , 0.999611699892566   , 0.546182784891012    , 0.0354392391275549   , 0                    ,
      1 , 1                  , 0.999804379743404   , 0.558391197313369    , 0.0393811953561564   , 0                    ,
      1 , 1                  , 0.999919683137743   , 0.570443625987811    , 0.043605133962085    , 0                    ,
      1 , 1                  , 0.999977477003604   , 0.582323875080213    , 0.0481207961547392   , 0                    ,
      1 , 1                  , 0.999997628269005   , 0.594015748756431    , 0.0529379231435176   , 0                    ,
      1 , 1                  , 1                   , 0.605503054739723    , 0.058066252860997    , 1.03373953304785e-09 ,
      1 , 1                  , 1                   , 0.61677370464036     , 0.0635117430162173   , 1.19668272694403e-06 ,
      1 , 1                  , 1                   , 0.627827618018029    , 0.0692692904062715   , 8.26991626437916e-06 ,
      1 , 1                  , 1                   , 0.638666899110818    , 0.0753317794500676   , 2.65385489447176e-05 ,
      1 , 1                  , 1                   , 0.649293652156815    , 0.0816920945665178   , 6.13203953608382e-05 ,
      1 , 1                  , 1                   , 0.659709981394105    , 0.0883431201745315   , 0.00011793327010562  ,
      1 , 1                  , 1                   , 0.669917991060777    , 0.0952777406930185   , 0.000201694987771943 ,
      1 , 1                  , 1                   , 0.679919785394917    , 0.10248884054089     , 0.000317923362952684 ,
      1 , 1                  , 1                   , 0.689717468634612    , 0.109969304137055    , 0.000471936210240724 ,
      1 , 1                  , 1                   , 0.699313145017949    , 0.117712015900426    , 0.000669051344228942 ,
      1 , 1                  , 1                   , 0.708708918783016    , 0.125709860249911    , 0.000914586579510216 ,
      1 , 1                  , 1                   , 0.717906894167899    , 0.133955721604421    , 0.00121385973067742  ,
      1 , 1                  , 1                   , 0.726909175410686    , 0.142442484382866    , 0.00157218861232345  ,
      1 , 1                  , 1                   , 0.735717866749463    , 0.151163033004157    , 0.00199489103904117  ,
      1 , 1                  , 1                   , 0.744335072422318    , 0.160110251887204    , 0.00248728482542346  ,
      1 , 1                  , 1                   , 0.752762896667338    , 0.169277025450916    , 0.0030546877860632   ,
      1 , 1                  , 1                   , 0.761003443722609    , 0.178656238114205    , 0.00370241773555327  ,
      1 , 1                  , 1                   , 0.769058817826219    , 0.18824077429598     , 0.00443579248848655  ,
      1 , 1                  , 1                   , 0.776931123216255    , 0.198023518415153    , 0.00526012985945592  ,
      1 , 1                  , 1                   , 0.784622464130804    , 0.207997354890632    , 0.00618074766305426  ,
      1 , 1                  , 1                   , 0.792134944807953    , 0.218155168141328    , 0.00720296371387444  ,
      1 , 1                  , 1                   , 0.799470669485789    , 0.228489842586152    , 0.00833209582650935  ,
      1 , 1                  , 1                   , 0.806631742402399    , 0.238994262644013    , 0.00957346181555184  ,
      1 , 1                  , 1                   , 0.81362026779587     , 0.249661312733822    , 0.0109323794955948   ,
      1 , 1                  , 1                   , 0.820438349904289    , 0.26048387727449     , 0.0124141666812312   ,
      1 , 1                  , 1                   , 0.827088092965744    , 0.271454840684927    , 0.0140241411870538   ,
      1 , 1                  , 1                   , 0.833571601218321    , 0.282567087384042    , 0.0157676208276555   ,
      1 , 1                  , 1                   , 0.839890978900107    , 0.293813501790746    , 0.0176499234176292   ,
      1 , 1                  , 1                   , 0.84604833024919     , 0.305186968323949    , 0.0196763667715678   ,
      1 , 1                  , 1                   , 0.852045759503656    , 0.316680371402562    , 0.0218522687040642   ,
      1 , 1                  , 1                   , 0.857885370901593    , 0.328286595445495    , 0.0241829470297112   ,
      1 , 1                  , 1                   , 0.863569268681088    , 0.339998524871657    , 0.0266737195631017   ,
      1 , 1                  , 1                   , 0.869099557080227    , 0.351809044099959    , 0.0293299041188286   ,
      1 , 1                  , 1                   , 0.874478340337098    , 0.363711037549313    , 0.0321568185114848   ,
      1 , 1                  , 1                   , 0.879707722689788    , 0.375697389638626    , 0.0351597805556631   ,
      1 , 1                  , 1                   , 0.884789808376383    , 0.387760984786811    , 0.0383441080659565   ,
      1 , 1                  , 1                   , 0.889726701634972    , 0.399894707412776    , 0.0417151188569578   ,
      1 , 1                  , 1                   , 0.89452050670364     , 0.412091441935434    , 0.0452781307432598   ,
      1 , 1                  , 1                   , 0.899173327820476    , 0.424344072773693    , 0.0490384615394556   ,
      1 , 1                  , 1                   , 0.903687269223565    , 0.436645484346464    , 0.0530014290601378   ,
      1 , 1                  , 1                   , 0.908064435150996    , 0.448988561072657    , 0.0571723511198996   ,
      1 , 1                  , 1                   , 0.912306929840855    , 0.461366187371183    , 0.0615565455333336   ,
      1 , 1                  , 1                   , 0.916416857531229    , 0.473771247660951    , 0.0661593301150328   ,
      1 , 1                  , 1                   , 0.920396322460206    , 0.486196626360871    , 0.07098602267959     ,
      1 , 1                  , 1                   , 0.924247428865872    , 0.498635207889855    , 0.0760419410415983   ,
      1 , 1                  , 1                   , 0.927972280986315    , 0.511079876666813    , 0.0813324030156503   ,
      1 , 1                  , 1                   , 0.931572983059621    , 0.523523517110655    , 0.086862726416339    ,
      1 , 1                  , 1                   , 0.935051639323877    , 0.53595901364029     , 0.0926382290582572   ,
      1 , 1                  , 1                   , 0.938410354017171    , 0.548379250674629    , 0.098664228755998    ,
      1 , 1                  , 1                   , 0.94165123137759     , 0.560777112632583    , 0.104946043324154    ,
      1 , 1                  , 1                   , 0.944776375643221    , 0.573145483933061    , 0.111488990577318    ,
      1 , 1                  , 1                   , 0.94778789105215     , 0.585477248994975    , 0.118298388330084    ,
      1 , 1                  , 1                   , 0.950687881842465    , 0.597765292237232    , 0.125379554397043    ,
      1 , 1                  , 1                   , 0.953478452252254    , 0.610002498078746    , 0.132737806592789    ,
      1 , 1                  , 1                   , 0.956161706519602    , 0.622181750938425    , 0.140378462731915    ,
      1 , 1                  , 1                   , 0.958739748882597    , 0.634295935235181    , 0.148306840629014    ,
      1 , 1                  , 1                   , 0.961214683579326    , 0.646337935387923    , 0.156528258098678    ,
      1 , 1                  , 1                   , 0.963588614847877    , 0.65830063581556     , 0.165048032955501    ,
      1 , 1                  , 1                   , 0.965863646926336    , 0.670176920937004    , 0.173871483014074    ,
      1 , 1                  , 1                   , 0.96804188405279     , 0.681959675171166    , 0.183003926088992    ,
      1 , 1                  , 1                   , 0.970125430465326    , 0.693641782936953    , 0.192450679994847    ,
      1 , 1                  , 1                   , 0.972116390402032    , 0.705216128653279    , 0.202217062546232    ,
      1 , 1                  , 1                   , 0.974016868100994    , 0.716675596739051    , 0.212308391557739    ,
      1 , 1                  , 1                   , 0.9758289678003      , 0.728013071613183    , 0.222729984843962    ,
      1 , 1                  , 1                   , 0.977554793738036    , 0.739221437694583    , 0.233487160219494    ,
      1 , 1                  , 1                   , 0.97919645015229     , 0.750293579402159    , 0.244585235498927    ,
      1 , 1                  , 1                   , 0.980756041281149    , 0.761222381154826    , 0.256029528496854    ,
      1 , 1                  , 1                   , 0.982235671362699    , 0.77200072737149     , 0.267825357027869    ,
      1 , 1                  , 1                   , 0.983637444635028    , 0.782621502471063    , 0.279978038906563    ,
      1 , 1                  , 1                   , 0.984963465336223    , 0.793077590872458    , 0.292492891947531    ,
      1 , 1                  , 1                   , 0.986215837704371    , 0.803361876994581    , 0.305375233965364    ,
      1 , 1                  , 1                   , 0.987396665977559    , 0.813467245256343    , 0.318630382774656    ,
      1 , 1                  , 1                   , 0.988508054393875    , 0.823386580076656    , 0.33226365619        ,
      1 , 1                  , 1                   , 0.989552107191404    , 0.833112765874429    , 0.346280372025987    ,
      1 , 1                  , 1                   , 0.990530928608234    , 0.842638687068571    , 0.360685848097212    ,
      1 , 1                  , 1                   , 0.991446622882453    , 0.851957228077996    , 0.375485402218268    ,
      1 , 1                  , 1                   , 0.992301294252146    , 0.861061273321611    , 0.390684352203746    ,
      1 , 1                  , 1                   , 0.993097046955403    , 0.869943707218327    , 0.406288015868241    ,
      1 , 1                  , 1                   , 0.993835985230308    , 0.878597414187055    , 0.422301711026344    ,
      1 , 1                  , 1                   , 0.99452021331495     , 0.887015278646705    , 0.438730755492649    ,
      1 , 1                  , 1                   , 0.995151835447416    , 0.895190185016187    , 0.455580467081749    ,
      1 , 1                  , 1                   , 0.995732955865792    , 0.903115017714412    , 0.472856163608236    ,
      1 , 1                  , 1                   , 0.996265678808166    , 0.910782661160289    , 0.490563162886703    ,
      1 , 1                  , 1                   , 0.996752108512624    , 0.918185999772729    , 0.508706782731744    ,
      1 , 1                  , 1                   , 0.997194349217254    , 0.92531791797064     , 0.527292340957951    ,
      1 , 1                  , 1                   , 0.997594505160144    , 0.932171300172935    , 0.546325155379916    ,
      1 , 1                  , 1                   , 0.997954680579379    , 0.938739030798525    , 0.565810543812234    ,
      1 , 1                  , 1                   , 0.998276979713046    , 0.945013994266318    , 0.585753824069496    ,
      1 , 1                  , 1                   , 0.998563506799234    , 0.950989074995224    , 0.606160313966296    ,
      1 , 1                  , 1                   , 0.998816366076029    , 0.956657157404157    , 0.627035331317226    ,
      1 , 1                  , 1                   , 0.999037661781518    , 0.962011125912021    , 0.64838419393688     ,
      1 , 1                  , 1                   , 0.999229498153788    , 0.967043864937732    , 0.67021221963985     ,
      1 , 1                  , 1                   , 0.999393979430927    , 0.971748258900198    , 0.692524726240729    ,
      1 , 1                  , 1                   , 0.999533209851021    , 0.976117192218328    , 0.715327031554111    ,
      1 , 1                  , 1                   , 0.999649293652157    , 0.980143549311034    , 0.738624453394587    ,
      1 , 1                  , 1                   , 0.999744335072422    , 0.983820214597226    , 0.762422309576751    ,
      1 , 1                  , 1                   , 0.999820438349904    , 0.987140072495814    , 0.786725917915195    ,
      1 , 1                  , 1                   , 0.99987970772269     , 0.990096007425708    , 0.811540596224513    ,
      1 , 1                  , 1                   , 0.999924247428866    , 0.992680903805818    , 0.836871662319297    ,
      1 , 1                  , 1                   , 0.99995616170652     , 0.994887646055054    , 0.862724434014141    ,
      1 , 1                  , 1                   , 0.999977554793738    , 0.99670911859233     , 0.889104229123637    ,
      1 , 1                  , 1                   , 0.999990530928608    , 0.99813820583655     , 0.916016365462378    ,
      1 , 1                  , 1                   , 0.999997194349217    , 0.99916779220663     , 0.943466160844957    ,
  };

  gsl_matrix_view expected = gsl_matrix_view_array(ix_data, sizeof(ix_data) / (6 * sizeof(double)), 6);

  munit_log_vector_dim(MUNIT_LOG_DEBUG, &taus.vector);
  munit_log_size(MUNIT_LOG_DEBUG, spline_df);
  munit_log_vector_dim(MUNIT_LOG_DEBUG, &knots.vector);
  munit_log_matrix_dim(MUNIT_LOG_DEBUG, &expected.matrix);

  gsl_matrix *result = gsl_matrix_alloc(taus.vector.size, (knots.vector.size - spline_df) + 1);
  munit_log_matrix_dim(MUNIT_LOG_DEBUG, result);

  munit_assert_size(expected.matrix.size1, ==, result->size1);
  munit_assert_size(expected.matrix.size2, ==, result->size2);

  qdm_ispline_matrix(result, &taus.vector, spline_df, &knots.vector);

  for (size_t i = 0; i < expected.matrix.size1; i++) {
    for (size_t j = 0; j < expected.matrix.size2; j++) {
      double e = gsl_matrix_get(&expected.matrix, i, j);
      double r = gsl_matrix_get(result, i, j);

      munit_logf(MUNIT_LOG_DEBUG, "(%zu, %zu) %0.17g ?= %0.17g (delta = %0.17g)", i, j, e, r, fabs(e - r));
      munit_assert_double_equal(e, r, 10);
    }
  }

  return MUNIT_OK;
}

MunitTest tests_ispline[] = {
  {
    "/qdm_ispline_matrix"   , // name
    test_qdm_ispline_matrix , // test
    NULL                    , // setup
    NULL                    , // tear_down
    MUNIT_TEST_OPTION_NONE  , // options
    NULL                    , // parameters
  },
  { NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL },
};

static const MunitSuite test_suite_ispline = {
  "/ispline"              , // name
  tests_ispline           , // tests
  NULL                    , // suites
  1                       , // iterations
  MUNIT_SUITE_OPTION_NONE , // options
};
