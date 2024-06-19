
    // std::vector<Vector3i> test;
    // Vector3i monomer(0,0,0);
      
    // // //Loop_move_test
    // monomer<<1,1,0;
    // test.push_back(monomer);
    // monomer<<1,0,0;
    // test.push_back(monomer);
    // monomer<<0,0,0;
    // test.push_back(monomer);
    // monomer<<0,1,0;
    // test.push_back(monomer);
    // monomer<<0,0,0;
    // test.push_back(monomer);
    // monomer<<0,0,0;
    // test.push_back(monomer);
    // monomer<<0,1,0;
    // test.push_back(monomer);
    // monomer<<0,0,0;
    // test.push_back(monomer);
    // monomer<<0,1,0;
    // test.push_back(monomer);
    // monomer<<0,0,0;
    // test.push_back(monomer);
    // monomer<<0,1,0;
    // test.push_back(monomer);
    // monomer<<0,0,0;
    // test.push_back(monomer);
    // monomer<<1,0,0;
    // test.push_back(monomer);
    // monomer<<2,0,0;
    // test.push_back(monomer); 
    // monomer<<3,0,0;
    // test.push_back(monomer);
    // monomer<<4,0,0;
    // test.push_back(monomer); 
    // monomer<<4,1,0;
    // test.push_back(monomer);  
    // for (int i=0;i<test.size();i++){
    //     for (int j=0;j<3;j++){
    //         printf("%d\t",test[i][j]);
    //     }      
    //     printf("\n") ;
    // }   
    // printf("change");
    // printf("\n") ;
    // loop_reduce_move(test, 9,0, 1)  ;
    // for (int i=0;i<test.size();i++){
    //     for (int j=0;j<3;j++){
    //         printf("%d\t",test[i][j]);
    //     }      
    //     printf("\n") ;
    // }    

    // //test
    // std::vector<int> tmp;
    // tmp.push_back(0);
    // tmp.push_back(1);
    // tmp.push_back(2);
    // tmp.push_back(3);
    // tmp.push_back(4);
    // tmp.push_back(5);
    // tmp.push_back(6);
    // tmp.push_back(7);
    // tmp.erase(tmp.begin()+1, tmp.begin()+3); 
    // for (auto use:tmp){
    //     std::cout<<use<<"\n";
    // }

    //==========================================
    // std::ofstream initial;
    // initial.open("E:/leilei/semester_1.2/data/gm12878_tad/initial.txt"); 
    // for (int i=0;i<polymer[0].size();i++){
    //     for (int j=0;j<3;j++){
    //         initial << polymer[0][i][j] <<"\t";
    //     } 
    //     initial  << std::endl;     
    //     // printf("\n") ;
    // }     
    //===========================================
    //     std::ofstream burn_in;
    // burn_in.open("E:/leilei/semester_1.2/data/gm12878_tad/burn_in.txt");
    
    // for (int i=0;i<polymer[0].size();i++){
    //     for (int j=0;j<3;j++){
    //         //printf("%d\t",polymer[0][i][j]);
    //         burn_in << polymer[0][i][j] <<"\t";
    //     }    
    //     burn_in  << std::endl;
    //     //printf("\n") ;
    // }  
    //============================================
    // for (auto l = 0; l < number_of_threads; l++) {
    //     std::ofstream burn_in;
    //     burn_in.open("E:/leilei/semester_1.2/data/gm12878_tad/burn_in_"+std::to_string(l)+".txt");
        
    //     for (int i=0;i<polymer[0].size();i++){
    //         for (int j=0;j<3;j++){
    //             //printf("%d\t",polymer[0][i][j]);
    //             burn_in << polymer[l][i][j] <<"\t";
    //         }    
    //         burn_in  << std::endl;
    //         //printf("\n") ;
    //     }  
    //     burn_in.close();
    // }

    //===============================================
        //     //print  normalize result
    //     std::ofstream burn_in;
    //     burn_in.open("E:/leilei/semester_1.2/data/gm12878_tad/normalize.txt");
        
    //     for (int i=0;i<pol_length;i++){
    //         for (int j=0;j<pol_length;j++){
    //             //printf("%d\t",polymer[0][i][j]);
    //             burn_in << tad_final_contacts[i][j] <<"\t";
    //         }    
    //         burn_in  << std::endl;
    //         //printf("\n") ;
    //     }  
    //     burn_in.close();  

    
    
    //test 把location——tad打印
    // for (auto elem :one_tad_locations[0]){
    //     for (auto use:elem.second){
    //         printf("%d\t",use);
    //     }
    //     printf("\n");
    // }
    // initial0.close();

    //    //print to test
    // std::ofstream contact;
    // contact.open("E:/leilei/semester_1.2/data/gm12878_tad/tmp/tad_contact.txt");
    
    // for (auto tmp:reference_one_tad_contact){

    //         //printf("%d\t",polymer[0][i][j]);
    //         contact << tmp.first.first<< "\t"<<tmp.first.second<<"\t"<<tmp.second <<"\n";   
    // }  
    // contact.close();



//     void loop_reduce_move(std::vector<Vector3i> &polymer, int site, int thread_num,int m,int update){
//     int flag=1;
//     int energy_change=0;
//     Vector3i prop_move1;
//     if (site>=mid_pol_length){flag=1;}
//     else{flag=-1;}
//     // printf("%d\n",flag);
//     if (site==pol_length-1 or site==0){
//         boundary_move(polymer,site,thread_num,m);
//     }
//     else if (((site+3)<pol_length) &&( site>3) && (polymer[(site+2*flag)] == polymer[site])){

//         reduce_count+=1;

//         Vector3i move_vector;
//         Vector3i pre_site;
//         int move_index=0;
//         int move_num=0;
//         for (int i=0;i<3;++i){
//             if ((polymer[site+3*flag][i]-polymer[site+2*flag][i])!=0){
//                 move_index=i;
//                 move_num=(polymer[site+3*flag][i]-polymer[site+2*flag][i]);
//                 break;                
//             }
//         }
//         move_vector=polymer[site+3*flag]-polymer[site+2*flag];
//         prop_move1=polymer[site]+ move_vector;
//         move_num=2*move_num;
//         move_vector=move_vector*2;
//         //burn in 阶段
//         if (update==0){
//             if(!(prop_move1[0]==polymer[site+flag][0] && prop_move1[1]==polymer[site+flag][1] && prop_move1[2]==polymer[site+flag][2])){
//                 update_contact_location(polymer, site+flag , pol_length, prop_move1,polymer[site+flag],thread_num,m);
//                 polymer[site+flag]=prop_move1;
//             }   

//             if (flag==1){
//                 for (int i=site+2;i<polymer.size();++i){
//                     pre_site=polymer[i];
//                     polymer[i][move_index]=polymer[i][move_index]+move_num;//prop_move1
                    
//                     update_contact_location(polymer, i , pol_length, polymer[i],pre_site,thread_num,m);
//                     // prop_move1=polymer[i]+move_vector;
//                     // polymer[i]=prop_move1;
//                 }            
//             }
//             else{
//                 for (int i=site-2;i>=0;--i){
//                     pre_site=polymer[i];
//                     polymer[i][move_index]=polymer[i][move_index]+move_num;
//                     update_contact_location(polymer, i , pol_length, polymer[i],pre_site,thread_num,m);
//                     // prop_move1=polymer[i]+move_vector;
//                     // polymer[i]=prop_move1;

//                 }             
//             }
//         }
//         //run stage 
//         else {   //queue 记录下prop move
//             std::queue<Vector3i> prop_moves;
//             prop_moves.push(prop_move1);
//             energy_change+=delta_E_other(polymer, site+1*flag, pol_length, prop_move1,thread_num);          
//             if (flag==1){
//                 for (int i=site+2;i<polymer.size();++i){
//                     // prop_move1=polymer[i]+move_vector;
                    
//                     prop_move1=polymer[i];
//                     prop_move1[move_index]=prop_move1[move_index]+move_num;
//                     prop_moves.push(prop_move1);
//                     //更新energy
//                     energy_change+=delta_E_other(polymer, i, pol_length, prop_move1,thread_num);
//                 }            
//             }
//             else{
//                 for (int i=site-2;i>=0;--i){
//                     // prop_move1=polymer[i]+move_vector;

//                     prop_move1=polymer[i];
//                     prop_move1[move_index]=prop_move1[move_index]+move_num;                    
//                     prop_moves.push(prop_move1);
//                     energy_change+=delta_E_other(polymer, i, pol_length, prop_move1,thread_num);
//                 }             
//             }
//             // if (energy_change!=0){
//             //     std::cout<<"energy_change"<<energy_change<<"\n";
//             // }
//             if (accept_move(energy_change)==1){

//                 prop_move1=prop_moves.front();
//                 if(!(prop_move1[0]==polymer[site+flag][0] && prop_move1[1]==polymer[site+flag][1] && prop_move1[2]==polymer[site+flag][2])){                    
//                     update_contact_location(polymer, site+flag , pol_length, prop_move1,polymer[site+flag],thread_num,m);
//                     polymer[site+flag]=prop_move1;               
//                 }   

//                 prop_moves.pop();
//                 if (flag==1){
//                     for (int i=site+2;i<polymer.size();++i){  
//                         prop_move1=prop_moves.front();
//                         update_contact_location(polymer, i , pol_length, prop_move1,polymer[i],thread_num,m);
//                         polymer[i]=prop_move1;
//                         prop_moves.pop();
//                     } 
//                 }   
//                 else{
//                     for (int i=site-2;i>=0;--i){
//                         prop_move1=prop_moves.front();
//                         update_contact_location(polymer, i , pol_length, prop_move1,polymer[i],thread_num,m);
//                         polymer[i]=prop_move1;
//                         prop_moves.pop();
//                     }
//                 }         
//             }
//         }        
//     }
// }



// void loop_add_move(std::vector<Vector3i> &polymer, int site, int thread_num,int m,int update){
//     int flag=1;
//     if (site>=mid_pol_length){flag=1;}
//     else{flag=-1;}
//     if (site==pol_length-1 or site==0){
//         boundary_move(polymer,site,thread_num,m);
//     }
//     else if (((site+2)<pol_length) && (site>=2) && ( polymer[site+2*flag] == 2*polymer[site+flag]-polymer[site]))//无loop,在同一条直线
//     {

//         add_count+=1;

//         int energy_change=0;
//         Vector3i move_vector;
//         Vector3i prop_move1;
//         Vector3i pre_site;
//         //确定loop的方向
//         int direction = unipath(gen);
//         int move_index=0;
//         int move_num=0;
//         Vector3i rotate_vector;
//         // move_vector=polymer[site]-polymer[site+2];//反向方向
//         for (int i=0;i<3;++i){
//             if ((polymer[site+2*flag][i]-polymer[site][i])!=0){
//                 move_index=i;
//                 move_num=(polymer[site+2*flag][i]-polymer[site][i]);
//                 break;                
//             }
//         }
//         for (int i=0;i<3;i++){
//             if (i==direction/3){rotate_vector[i]=(-2*(direction%2)+1);}//
//             else{rotate_vector[i]=0;}
//         }

//         move_vector=polymer[site+2*flag]-polymer[site]; 
//         prop_move1=polymer[site]+rotate_vector;
//         if (update==0){

//             if(!(prop_move1[0]==polymer[site+flag][0] && prop_move1[1]==polymer[site+flag][1] && prop_move1[2]==polymer[site+flag][2])){
//                 update_contact_location(polymer, site+flag , pol_length, prop_move1,polymer[site+flag],thread_num,m);
//                 polymer[site+flag]=prop_move1;
//             }   
    
//             if (flag==1){
//                 for (int i=site+2;i<polymer.size();++i){
//                     pre_site=polymer[i];
//                     polymer[i][move_index]=polymer[i][move_index]-move_num;//prop_move1
//                     update_contact_location(polymer, i , pol_length, polymer[i],pre_site,thread_num,m);
//                     // prop_move1=polymer[i]-move_vector;
//                     // update_contact_location(polymer, i , pol_length, prop_move1,polymer[i],thread_num,m);
//                     // polymer[i]=prop_move1;
                    
//                 }            
//             }
//             else{
//                 for (int i=site-2;i>=0;--i){
//                     pre_site=polymer[i];
//                     polymer[i][move_index]=polymer[i][move_index]-move_num;//prop_move1
//                     update_contact_location(polymer, i , pol_length, polymer[i],pre_site,thread_num,m);
//                     // prop_move1=polymer[i]-move_vector;
//                     // update_contact_location(polymer, i , pol_length, prop_move1,thread_num,m);
//                 }             
//             }    
//         }
//         else {   //queue 记录下prop move
//             std::queue<Vector3i> prop_moves;
//             prop_moves.push(prop_move1);
//             energy_change+=delta_E_other(polymer, site+flag, pol_length, prop_move1,thread_num);       
//             if (flag==1){
//                 for (int i=site+2;i<polymer.size();++i){
//                     // prop_move1=polymer[i]-move_vector;
//                     prop_move1=polymer[i];
//                     prop_move1[move_index]=prop_move1[move_index]+move_num;
//                     prop_moves.push(prop_move1);
//                     //更新energy
//                     energy_change+=delta_E_other(polymer, i, pol_length, prop_move1,thread_num);
//                 }            
//             }
//             else{
//                 for (int i=site-2;i>=0;--i){
//                     // prop_move1=polymer[i]-move_vector;
//                     prop_move1=polymer[i];
//                     prop_move1[move_index]=prop_move1[move_index]+move_num;
//                     prop_moves.push(prop_move1);
//                     //更新energy
//                     energy_change+=delta_E_other(polymer, i, pol_length, prop_move1,thread_num);                    
//                 }             
//             }
//             if (accept_move(energy_change)==1){

//                 prop_move1=prop_moves.front();
//                 if(!(prop_move1[0]==polymer[site+flag][0] && prop_move1[1]==polymer[site+flag][1] && prop_move1[2]==polymer[site+flag][2])){
//                     update_contact_location(polymer, site+flag , pol_length, prop_move1,polymer[site+flag],thread_num,m);
//                     polymer[site+flag]=prop_move1;
//                 }                  
//                 prop_moves.pop();

//                 if (flag==1){
//                     for (int i=site+2;i<polymer.size();++i){  
//                         prop_move1=prop_moves.front();
//                         update_contact_location(polymer, i , pol_length, prop_move1,polymer[i],thread_num,m);
//                         polymer[i]=prop_move1;
//                         prop_moves.pop();
//                     } 
//                 }   
//                 else{
//                     for (int i=site-2;i>=0;--i){
//                         prop_move1=prop_moves.front();
//                         update_contact_location(polymer, i , pol_length, prop_move1,polymer[i],thread_num,m);
//                         polymer[i]=prop_move1;
//                         prop_moves.pop();
//                     }
//                 }         
//             }
//         }
//     }
// }
