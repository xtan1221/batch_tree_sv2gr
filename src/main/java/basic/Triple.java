/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package basic;

/**
 *
 * @author tanxu
 */
public class Triple<L,M,R>{
        L l;
        M m;
        R r;
        public Triple(L l,M m, R r){
            this.l = l;
            this.m = m;
            this.r = r;
        }
        
        public L getLeft(){
            return l;
        }
        public M getMiddle(){
            return m;
        }
        public R getRight(){
            return r;
        }
        
        
        
        @Override
        public String toString(){
            StringBuilder sb = new StringBuilder();
            sb.append(l).append(",").append(m).append(",").append(r);
            return sb.toString();
        }
    }
