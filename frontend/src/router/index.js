import { createRouter, createWebHistory } from 'vue-router'

import Layout from '@/views/Layout/index.vue'

// Views are code-split so the first paint only downloads the shell and the
// requested page; Element Plus components are auto-imported per chunk.
const router = createRouter({
  history: createWebHistory(import.meta.env.BASE_URL),
  scrollBehavior(to, from, savedPosition) {
    if (savedPosition) return savedPosition
    if (to.path === from.path) return false
    return { top: 0 }
  },
  routes: [
    {
      path: '/',
      component: Layout,
      redirect: '/Home',
      children: [
        {
          path: 'Home',
          name: 'Home',
          component: () => import('@/views/Home/index.vue'),
        },
        {
          path: 'Dataset',
          name: 'Dataset',
          component: () => import('@/views/Dataset/index.vue'),
        },
        {
          path: 'Gene',
          name: 'Gene',
          component: () => import('@/views/Gene/index.vue'),
        },
        {
          path: 'KEGG',
          name: 'KEGG',
          component: () => import('@/views/KEGG/index.vue'),
        },
        {
          path: 'Information',
          name: 'Information',
          component: () => import('@/views/Documentation/index.vue'),
        },
        {
          path: 'TeamInformation',
          name: 'TeamInformation',
          component: () => import('@/views/TeamInformation/index.vue'),
        },
        {
          path: ':pathMatch(.*)*',
          name: 'NotFound',
          component: () => import('@/views/NotFound/index.vue'),
        },
      ],
    },
    {
      path: '/login',
      component: () => import('@/views/Login/index.vue'),
    },
  ],
})

export default router
