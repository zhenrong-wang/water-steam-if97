# Water-and-Steam-Property-Calculation
High-resolution water and steam property calculation, based on IAPWS-IF 97

Developer: Zhenrong WANG. zhenrongwang@live.com 
WeChat: K495458966


1.	本程序名为_WSPC.exe（Water and Steam Properties Calculation)，版本为1.1beta。技术参照：IAPWS IF97。http://www.iapws.org/

2.	程序功能是根据给定压力（p）、温度（t）、密度（r）、内能（u）、比焓（h）、比熵（s）六个参数以及干度（x）中的两个，计算出一系列物性参数。

3.	程序输入文件：_input.dat。格式为每行三个数据，第一个为整数，表示计算类型（给定参数）：1-给定pt，2-给定pr，3-给定pu，4-给定ph，5-给定ps，6-给定tr，7-给定tu，8-给定th，9-给定ts，10-给定hs，11-给定px，12-给定tx（x为干度）。第二、三个参数为浮点数（可以使用C语言认可的科学计数法，如103=1.03e2。请见范例文件）。请务必注意参数顺序，严格按照ptruhs的先后顺序。各参数之间用逗号隔开。所有行之后不得再有任何空格、回车等。

4.	各个参数的单位需如下：p-Pa, t-K, r-kg/m3, u-J/kg, h-J/kg, s-J/(kg*K)，请务必注意。

5.	程序输出文件：_properties.dat，格式见范例文件。

